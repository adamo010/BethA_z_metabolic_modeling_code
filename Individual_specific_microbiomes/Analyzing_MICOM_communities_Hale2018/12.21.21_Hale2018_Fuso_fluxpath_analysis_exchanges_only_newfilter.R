library("ggplot2")
library("readxl")
library("dplyr")
library("tidyverse")
library("ggpubr")
library("ape")
library("lme4")
library("gplots")
library("plotly")
library("tidyr")
library("vegan")
library("data.table")
library("stringr")

#the goal here is to calculate p-values for the significant differences between flux pathways to use for correlating with gene expression
#differences. Also, we want to abundance filter flux pathways that do not appear in at least half the number of samples 
#(before running stats on identifying fluxes with significantly different activity between tumor and normal samples).

#this is different from previous kinds of analyses because we collapsed differently here; only looking at exchange fluxes (named EX_), and
#did not collapse into a single value for each community (i.e. taxon-specific values are preserved). 

#starting point is 10.26.21.21_Hale2018_fluxpath_analysis_exchanges_only.R
#and 11.08.21_Burns2015_Fuso_fluxpath_analysis_exchanges_only.R
#and 11.09.21_Hale2018_Fuso_fluxpath_analysis_exchanges_only.R- here, we're using different (cleaner) filtering

#first, we'll use the original OTU table to create a list of OTU_ids to keep; these are the fuso OTU ids

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Hale2018/") 

##########Step 0: finding fuso only otu_ids (NEW/moded)
#pull out fuso only otu_ids (there are multiple for this dataset)
taxonomy <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Hale2018/Hale_2018_OTUs_with_matched_taxonomy_and_proxy_OTUids.csv")
fuso_only_samples <- filter(taxonomy, Genus == "Fusobacterium")
fuso_otu_id <- fuso_only_samples$OTU_ID
#clean up
rm(taxonomy, fuso_only_samples)

##########Step 1: import taxonomy data and use to get a list of Fuso-specific OTU ids
#use a post-filtered abundance table:
filtered_tax <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Hale2018/09.01.21_Hale2018_indiv_spp_GRs_collapsed_to_genus.csv")
fuso_only_tax <- filter(filtered_tax, Genus == "Fusobacterium") #filter to only include fuso samples
#set up a sorting df
fuso_pos_samples <- subset(fuso_only_tax, select = c(host_subject_id, SampleID)) #dump unnecessary cols for sorting
paired_fuso_PBIDs <- fuso_pos_samples[duplicated(fuso_pos_samples[,"host_subject_id"]) | duplicated(fuso_pos_samples[,"host_subject_id"], fromLast=TRUE),] #remove unpaired PBIDs
paired_fuso_PBIDs$host_subject_id <- as.factor(paired_fuso_PBIDs$host_subject_id)
paired_fuso_PBIDs2 <- unique(paired_fuso_PBIDs$host_subject_id) #a list of paired Fuso+ PBIDs
#clean up
rm(filtered_tax, fuso_pos_samples, paired_fuso_PBIDs, fuso_only_tax)

##########Step 2: import data and adjust factors
analyzed_fluxpaths <- read_csv("10.13.21_Hale2018_fluxes_with_metadata_longform.csv")
analyzed_fluxpaths$host_subject_id <- as.factor(analyzed_fluxpaths$host_subject_id)
analyzed_fluxpaths$normal_adjacent_or_tumor_tissue_specimen <- as.factor(analyzed_fluxpaths$normal_adjacent_or_tumor_tissue_specimen)
#drop unnecessary columns
analyzed_fluxpaths <- subset(analyzed_fluxpaths, select = -c(...1, fluxpath_formula, fluxpath_subsystem, sampleIDs_from_read_counts))
#edit column names to match Burns because I am lazy and don't want to think too much about rerunning this code again.
analyzed_fluxpaths <- analyzed_fluxpaths %>% rename(Patient_Blind_ID = host_subject_id, Description = normal_adjacent_or_tumor_tissue_specimen)
unique_fluxpaths <- unique(analyzed_fluxpaths$fluxpath_name) #324 fluxpaths
#clean up
rm(unique_fluxpaths)

#########NEW Step 3: filtering by fuso only samples
fuso_fluxpaths_only <- analyzed_fluxpaths[analyzed_fluxpaths$OTU_ID %in% fuso_otu_id, ] #filter to only include fuso OTU_ids
#use the sorting DF created in step 1
fuso_fluxpaths_only_filtered <- fuso_fluxpaths_only[fuso_fluxpaths_only$Patient_Blind_ID %in% paired_fuso_PBIDs2, ] #filter by new list
fuso_PBIDs_unique2 <- unique(fuso_fluxpaths_only_filtered$Patient_Blind_ID)
#clean up
rm(fuso_otu_id, analyzed_fluxpaths, fuso_fluxpaths_only, paired_fuso_PBIDs2)
#wonderful. okay, the key here was (for some dumb fuck reason) to NOT set Patient_Blind_ID as factor before changing.
fuso_fluxpaths_only_filtered$Patient_Blind_ID <- as.factor(fuso_fluxpaths_only_filtered$Patient_Blind_ID)

#how many fluxpaths do we have? 274
fuso_fluxpaths_unique <- unique(fuso_fluxpaths_only_filtered$fluxpath_name)
rm(fuso_fluxpaths_unique)

#########Step 3a: adding across OTU_ids
#because we have multiple fuso_ids in this dataset, we need to collapse these. I don't want to fuck around with multiple OTUs per sample
#here, we would like to sum all flux_values across OTU_ids, within fluxpath_name/sample_ID combos
#first, drop unnecessary columns (i.e. what we're adding across)
sum_fuso_fluxpaths <- subset(fuso_fluxpaths_only_filtered, select = -c(OTU_ID))
#create a new sorting column to get unique values for each Sample_ID/Fluxpath_name combo
sum_fuso_fluxpaths <- sum_fuso_fluxpaths %>% 
  mutate(sorting_col = paste0(Sample_ID, "_", fluxpath_name)) 
#Drop NAN values; converting to zeroes screws with averages
#sum_fuso_fluxpaths_noNANs <- sum_fuso_fluxpaths[!is.na(sum_fuso_fluxpaths$flux_value),]
#okay, here is another source of contention. I don't want to lose fluxpaths. And really, isn't zero basically NAN, since we've already
#filtered by Fuso-specific samples?
#NEW 12.21.21: replace NANs wtih zeroes
sum_fuso_fluxpaths$flux_value <- as.numeric(sum_fuso_fluxpaths$flux_value)
sum_fuso_fluxpaths$flux_value[is.na(sum_fuso_fluxpaths$flux_value)] <- 0

#then sum across sorting col:
sum_fuso_fluxpaths_noNANs_collapsed <- sum_fuso_fluxpaths %>%        # Specify data frame
  group_by(sorting_col, fluxpath_name, Sample_ID, fluxpath_description, Patient_Blind_ID, Description) %>%    # Specify group indicator (columns to keep)
  summarise_at(vars(flux_value),                   # Specify column
               list(sum_fluxpath_amount = sum)) %>% # Specify function
  ungroup() #do this so select (dropping columns) can happen later
#nice. End up with a dataframe where each row is a unique subsystem/sample average of fluxes from pathways within that subsystem
#how many unique fluxes are we including here? 
sum_fluxpaths_unique <- unique(sum_fuso_fluxpaths_noNANs_collapsed$fluxpath_name)

#rename for latter consistency in code
sum_fuso_fluxpaths_noNANs_collapsed <- sum_fuso_fluxpaths_noNANs_collapsed %>% rename(flux_value = sum_fluxpath_amount)

#clean up
rm(fuso_fluxpaths_only_filtered, sum_fuso_fluxpaths, sum_fluxpaths_unique)

#########Step 4: collapsing by trace minerals: borrowed from 09.20.21_Burns2015_inputoutput_analysis_plus_traceminerals.R
#there are several metabolites (EX_ca2(e), EX_cobalt2(e), EX_cu2(e), EX_mg2(e), EX_mn2(e), EX_zn2(e), EX_k(e), EX_cl(e), EX_so4(e))
#which have identical fluxes in all communities. The goal is to collapse these into a single group, called "Trace minerals"
#first, create a vector of these metabolites
trace_minerals <- c("EX_ca2(e)", "EX_cobalt2(e)", "EX_cu2(e)", "EX_mg2(e)", "EX_mn2(e)", "EX_zn2(e)", "EX_k(e)", "EX_cl(e)", "EX_so4(e)")
#Figure out some way to count the number of trace minerals with nonzero values within each sample. 
#As long as it's always 10 or 0, can collapse. 
#first, create a dataframe that ONLY contains trace minerals
fuso_fluxpaths_trace_minerals_only <- sum_fuso_fluxpaths_noNANs_collapsed[sum_fuso_fluxpaths_noNANs_collapsed$fluxpath_name %in% trace_minerals, ]

#use dplyr to round to 6 significant digits.
fuso_fluxpaths_trace_minerals_only <- fuso_fluxpaths_trace_minerals_only %>% 
  mutate(flux_value=signif(flux_value, 6)) 

#then, create a filtering_DF that counts the number of unique values within each sample_ID_ish
trace_minerals_filtering_df <- fuso_fluxpaths_trace_minerals_only %>%  #create a new dataframe by filtering the old one
  group_by(Sample_ID) %>% #group by the column of interest
  summarize(n = n_distinct(flux_value)) #count the number of occurrences of each unique value in metabolite_amount and store as variable 'n'
#okay, now we have all 1s, with a few 3s. This is actually good. It means within each sample_ID_ish, there is only 1-3 unique
#metabolite_amounts. The differences are so small as to be ignore-able

#time to collapse! for each sample_ID_ish, only keep ONE of the rows where metabolite_names values are in the trace_minerals vector
#and rename metabolite_name to trace_mineral.
fuso_fluxpaths_collapsed <-data.frame(sum_fuso_fluxpaths_noNANs_collapsed) #copy the original dataframe
#first, rename "EX_cobalt2(e)" to "exchange of trace elements"
fuso_fluxpaths_collapsed$fluxpath_name[fuso_fluxpaths_collapsed$fluxpath_name== "EX_cobalt2(e)"] <- "Trace_elements_exchange"
fuso_fluxpaths_collapsed$fluxpath_description[fuso_fluxpaths_collapsed$fluxpath_description== "Co2+ exchange"] <- "Exchange of nine trace elements"
#then, remove all rows where trace_minerals vector values are in the metabolite_name column
fuso_fluxpaths_collapsed<- fuso_fluxpaths_collapsed[!(fuso_fluxpaths_collapsed$fluxpath_name %in% trace_minerals),]
#number of unique values
fuso_fluxpaths_collapsed_unique <- unique(fuso_fluxpaths_collapsed$fluxpath_name)

#clean up
rm(sum_fuso_fluxpaths_noNANs_collapsed, trace_minerals, fuso_fluxpaths_trace_minerals_only, trace_minerals_filtering_df, fuso_fluxpaths_collapsed_unique)

#########Step 5: remove samples where tumor and normal fluxes = 0
fuso_fluxpaths_wide <- select(fuso_fluxpaths_collapsed, -c(Sample_ID, sorting_col)) #drop unnecessary columns
#first, convert to wide form
fuso_fluxpaths_wide <- fuso_fluxpaths_wide %>% 
  spread(Description, flux_value)
#then, drop all rows where both tumor and normal values are nan.
fuso_fluxpaths_wide_noNANs <- fuso_fluxpaths_wide[!(is.na(fuso_fluxpaths_wide$tumor) & is.na(fuso_fluxpaths_wide$normal)),]
fuso_fluxpaths_noNANs_unique <- unique(fuso_fluxpaths_wide_noNANs$fluxpath_name) #count the number of metabolites
#replace remaining NA values with zeroes
fuso_fluxpaths_wide_noNANs$tumor[is.na(fuso_fluxpaths_wide_noNANs$tumor)] = 0
fuso_fluxpaths_wide_noNANs$normal[is.na(fuso_fluxpaths_wide_noNANs$normal)] = 0
#create a difference column
fuso_fluxpaths_wide_noNANs$difference <- (fuso_fluxpaths_wide_noNANs$tumor - fuso_fluxpaths_wide_noNANs$normal)
#remove any rows where differences are 0
fuso_fluxpaths_wide_noNANs_zerofree <- fuso_fluxpaths_wide_noNANs[(fuso_fluxpaths_wide_noNANs$difference !=0),]
#NOTE: for this dataset, a few rows are removed but no fluxpaths
#how many unique fluxes are we including here? 
fuso_fluxpaths_zerofree_unique <- unique(fuso_fluxpaths_wide_noNANs_zerofree$fluxpath_name)

#cleanup: note to KEEP fuso_fluxpaths_collapsed
rm(fuso_fluxpaths_wide, fuso_fluxpaths_wide_noNANs, fuso_fluxpaths_noNANs_unique, fuso_fluxpaths_zerofree_unique)

########Step 6: filter by flux pathways that are active in ALL Fuso-containing samples
#first, find out how many samples have Fuso.
num_fuso_pos_samples <- unique(fuso_fluxpaths_collapsed$Sample_ID) #46 samples, here.

#first, create a filtering dataframe that contains counts of each fluxpath_name incidence in the dataframe
fuso_fluxpaths_filtering_df <- fuso_fluxpaths_wide_noNANs_zerofree %>%  #create a new dataframe by filtering the old one
  mutate(count_of_nonzeros = rowSums(.[4:5]!=0)) %>% #create a new column, count_of_nonzeroes, that counts # nonzero values in tumor and normal columns
  group_by(fluxpath_name) %>% #group by the column of interest
  summarize(n = sum(count_of_nonzeros)) %>% #sums the count_of_nonzeros values within each fluxpath_subsystem and saves the value as n
  mutate(n) #add n as a new column in the dataframe
#create a new dataframe that merges the old and new dataframes
fuso_fluxpaths_filter_added <- dplyr::inner_join(fuso_fluxpaths_wide_noNANs_zerofree, fuso_fluxpaths_filtering_df, by= "fluxpath_name")
#now, remove all rows where n (the name of the count column) is less than ... 46 is 100%, 23 is 50%
fuso_fluxpaths_abundance_filtered <- filter(fuso_fluxpaths_filter_added, n >= 23)
#count number of unique metabolites left
fuso_fluxpaths_filtered_unique <- unique(fuso_fluxpaths_abundance_filtered$fluxpath_name)
#cleanup
rm(fuso_fluxpaths_filtering_df,fuso_fluxpaths_filter_added, fuso_fluxpaths_only_filtered, num_fuso_pos_samples, fuso_fluxpaths_filtered_unique, fuso_fluxpaths_wide_noNANs_zerofree)

##########Step 7: split into multiple dataframes by metabolite and do stats on each metabolite
#first, edit the dataframe; drop a couple of columns
fuso_fluxpaths_paired_only <- select(fuso_fluxpaths_abundance_filtered, -c(difference, n, fluxpath_description))
#convert to long form
fuso_fluxpaths_paired_long <- reshape2::melt(data= fuso_fluxpaths_paired_only,
                                             id.vars= c("Patient_Blind_ID", "fluxpath_name"),
                                             variable.name = "Description",
                                             value.name = "fluxpath_activity")
#split up by different subsystems and see where we get. 
fuso_split_analyzed_fluxpaths_by_flux <- split(fuso_fluxpaths_paired_long, 
                                               with(fuso_fluxpaths_paired_long, interaction(fluxpath_name)), drop = TRUE)
fuso_by_flux_stats <- lapply(fuso_split_analyzed_fluxpaths_by_flux, function(df){wilcox.test(fluxpath_activity~Description, data=df, exact= FALSE, paired= TRUE)})
#make lists of p values
by_flux_pvals <- c() #initialize list
for (elem in fuso_by_flux_stats){
  new_value = elem$p.value
  by_flux_pvals <- c(by_flux_pvals, new_value)}
rm(new_value)
#make a list of q-values
by_flux_qvals <- p.adjust(by_flux_pvals, method = "fdr")
#make a list of genus_ids
flux_ids <- c()
for(elem in fuso_split_analyzed_fluxpaths_by_flux){
  new_value = elem$fluxpath_name[1]
  flux_ids <- c(flux_ids, new_value)}
#merge all lists together. 
fuso_fluxpath_statistics <- data.frame(flux_ids, by_flux_pvals, by_flux_qvals)
#clean up and remove extra variables
rm(elem, fuso_split_analyzed_fluxpaths_by_flux, by_flux_pvals, by_flux_qvals, flux_ids, new_value)

###########Step 8: sort and save
fuso_fluxpath_statistics <- fuso_fluxpath_statistics[order(fuso_fluxpath_statistics$by_flux_qvals),] #re-order by q-value
fuso_fluxpath_statistics <- fuso_fluxpath_statistics %>% rename("fluxpath_name"="flux_ids" )
#I would also like to add back the descriptions, please.
fluxpath_key <- subset(fuso_fluxpaths_collapsed, select = c(fluxpath_name, fluxpath_description) )
fluxpath_key <- fluxpath_key %>% distinct(fluxpath_name, .keep_all = TRUE)
#create a new dataframe that merges the old and new dataframes
fuso_fluxpath_statistics2 <- dplyr::inner_join(fuso_fluxpath_statistics, fluxpath_key, by= "fluxpath_name")
write.csv(as.data.frame(fuso_fluxpath_statistics2), file="12.21.21_Fusobacterium_exchange_fluxes_only_stats_Hale2018_data.csv")

##########Step 9: prepare data for graphing 
#add a new column called diff. tumor-normal
fuso_fluxpaths_paired_only$diff <- (fuso_fluxpaths_paired_only$tumor - fuso_fluxpaths_paired_only$normal)
#NOW we can add a new column with useful binning of differences for graphing purposes
fuso_fluxpaths_paired_wide <- fuso_fluxpaths_paired_only %>%
  mutate(Difference_direction = case_when(
    diff > 0 ~ "tumor > normal",
    diff < 0 ~ "tumor < normal",
    diff == 0 ~ "Zero change",))
#I would also like to add back the descriptions, please.
fluxpath_key <- subset(fuso_fluxpaths_collapsed, select = c(fluxpath_name, fluxpath_description) )
fluxpath_key <- fluxpath_key %>% distinct(fluxpath_name, .keep_all = TRUE)
#create a new dataframe that merges the old and new dataframes
fuso_fluxpaths_paried_wide2 <- dplyr::inner_join(fuso_fluxpaths_paired_wide, fluxpath_key, by= "fluxpath_name")
write.csv(as.data.frame(fuso_fluxpaths_paried_wide2), file="12.21.21_Fusobacterium_exchange_fluxes_only_alldata_Hale2018_data.csv")

#step 10: graphing!
#Acetate
ggpaired(subset(fuso_fluxpaths_paried_wide2, fluxpath_name %in% c("EX_ac(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Fusobacterium acetate exchange flux,\n Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.000744, q= 0.0360")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Sum of Fusobacterium fluxes,\nmmol/(gDW * h)", limits=c(0,5000)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Carbon dioxide
ggpaired(subset(fuso_fluxpaths_paried_wide2, fluxpath_name %in% c("EX_co2(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Fusobacterium carbon dioxide exchange flux,\n Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.00159, q= 0.0360")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Sum of Fusobacterium fluxes,\nmmol/(gDW * h)", limits=c(0,6000)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Glucose
ggpaired(subset(fuso_fluxpaths_paried_wide2, fluxpath_name %in% c("EX_glc_D(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Fusobacterium glucose exchange flux,\n Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.00128, q= 0.0360")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Sum of Fusobacterium fluxes,\nmmol/(gDW * h)", limits=c(-500,500)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

