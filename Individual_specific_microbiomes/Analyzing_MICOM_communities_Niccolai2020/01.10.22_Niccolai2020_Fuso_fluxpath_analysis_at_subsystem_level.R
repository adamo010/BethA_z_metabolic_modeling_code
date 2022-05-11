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

#this analysis is specifically focused on idenitying Fuso flux subsystem differences between tumor and normal samples
#I've already done ALL subsystems, and Fuso exchanges only. Now we're doing fuso subsystems. You know, out of desperation.

#based on combining 01.07.22_Hale2018_Fuso_fluxpath_analysis_at_subsystem_level.R 
#and 12.21.21_Niccolai2020_Fuso_fluxpath_analysis_exchanges_only_newfilter.R

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Niccolai2020/") 

##########Step 2: import data and adjust factors
analyzed_fluxpaths <- read_csv("01.10.22_Niccolai2020_Fuso_fluxes_with_metadata_longform.csv") 
analyzed_fluxpaths$Patient_Blind_ID <- as.factor(analyzed_fluxpaths$Patient_Blind_ID)
analyzed_fluxpaths <- analyzed_fluxpaths %>% rename(fluxpath_name= fluxpath_id, 
                                                    fluxpath_subsystem= subsystem,
                                                    flux_value = fluxpath_values,
                                                    fluxpath_description = description)
#how many unique fluxes are we including here? 
fluxpaths_unique <- unique(analyzed_fluxpaths$fluxpath_name) #7522 fluxpaths
#how many PBIDs?
PBIDs_unique <- unique(analyzed_fluxpaths$Patient_Blind_ID) #36 PBIDs
#how many subsystems?
subsystems_unique <- unique(analyzed_fluxpaths$fluxpath_subsystem) #105 subsystems
#clean up
rm(fluxpaths_unique, PBIDs_unique, subsystems_unique)

#########Step 3a: adding across OTU_ids
#because we have multiple fuso_ids in this dataset, we need to collapse these. I don't want to fuck around with multiple OTUs per sample
#here, we would like to sum all flux_values across OTU_ids, within fluxpath_name/sample_ID combos
#first, drop unnecessary columns (i.e. what we're adding across)
sum_fuso_fluxpaths <- subset(analyzed_fluxpaths, select = -c(OTU_ID))
#create a new sorting column to get unique values for each Sample_ID/Fluxpath_name combo
sum_fuso_fluxpaths <- sum_fuso_fluxpaths %>% 
  mutate(sorting_col = paste0(SampleID, "_", fluxpath_name)) 
#Drop NAN values; converting to zeroes screws with averages
#sum_fuso_fluxpaths_noNANs <- sum_fuso_fluxpaths[!is.na(sum_fuso_fluxpaths$flux_value),]
#okay, here is another source of contention. I don't want to lose fluxpaths. And really, isn't zero basically NAN, since we've already
#filtered by Fuso-specific samples?
#NEW 12.21.21: replace NANs wtih zeroes
sum_fuso_fluxpaths$flux_value <- as.numeric(sum_fuso_fluxpaths$flux_value)
sum_fuso_fluxpaths$flux_value[is.na(sum_fuso_fluxpaths$flux_value)] <- 0

#then sum across sorting col:
sum_fuso_fluxpaths_noNANs_collapsed <- sum_fuso_fluxpaths %>%        # Specify data frame
  group_by(sorting_col, fluxpath_name, SampleID, fluxpath_description, Patient_Blind_ID, Description) %>%    # Specify group indicator (columns to keep)
  summarise_at(vars(flux_value),                   # Specify column
               list(sum_fluxpath_amount = sum)) %>% # Specify function
  ungroup() #do this so select (dropping columns) can happen later
#nice. End up with a dataframe where each row is a unique subsystem/sample average of fluxes from pathways within that subsystem
#how many unique fluxes are we including here? 
sum_fluxpaths_unique <- unique(sum_fuso_fluxpaths_noNANs_collapsed$fluxpath_name)

#rename for latter consistency in code
sum_fuso_fluxpaths_noNANs_collapsed <- sum_fuso_fluxpaths_noNANs_collapsed %>% rename(flux_value = sum_fluxpath_amount)

#clean up
rm(sum_fuso_fluxpaths, sum_fluxpaths_unique)

#########Step 4: collapse this dataset by subsystem; average across fluxpaths within a subsystem within a sample.
#drop unnecessary columns
analyzed_fluxpaths_subset <- subset(analyzed_fluxpaths, select = -c(...1, abbreviation, fluxpath_description, formula))
#create a new sorting column to get unique values for each Sample_ID/Fluxpath_subsystem combo
analyzed_fluxpaths_subset <- analyzed_fluxpaths_subset %>% 
  mutate(sorting_col = paste0(SampleID, "_", fluxpath_subsystem)) 
#AAAA I think I have to replace NANs with 0s. Otherwise the averages won't work 
#analyzed_fluxpaths_subset$fluxpath_amount[is.na(analyzed_fluxpaths_subset$fluxpath_amount)] = 0
#hmmmm, is it better to convert to zeroes or drop??? Probably drop. Let's go with that.
analyzed_fluxpaths_subset_noNANs <- analyzed_fluxpaths_subset[!is.na(analyzed_fluxpaths_subset$flux_value),]
#then average across sorting col:
analyzed_fluxpaths_subset_averaged <- analyzed_fluxpaths_subset_noNANs %>%        # Specify data frame
  group_by(sorting_col, SampleID, Description, Patient_Blind_ID, fluxpath_subsystem) %>%    # Specify group indicator (columns to keep)
  summarise_at(vars(flux_value),                   # Specify column
               list(av_fluxpath_amount = mean)) %>% # Specify function
  ungroup() #do this so select (dropping columns) can happen later
#nice. End up with a dataframe where each row is a unique subsystem/sample average of fluxes from pathways within that subsystem
#how many unique fluxes are we including here? 
fluxpaths_averaged_unique <- unique(analyzed_fluxpaths_subset_averaged$fluxpath_subsystem) #85 subsystems
#clean up
rm(fluxpaths_averaged_unique, analyzed_fluxpaths_subset_noNANs, analyzed_fluxpaths_subset)

#########Step 3: filter out zeroes
analyzed_fluxpaths_wide <- select(analyzed_fluxpaths_subset_averaged, -c(sorting_col, SampleID)) #drop unnecessary columns
#first, convert to wide form
analyzed_fluxpaths_wide <- analyzed_fluxpaths_wide %>% 
  spread(Description, av_fluxpath_amount)
#then, drop all rows where both tumor and normal values are nan.
analyzed_fluxpaths_wide_noNANs <- analyzed_fluxpaths_wide[!(is.na(analyzed_fluxpaths_wide$tumor) & is.na(analyzed_fluxpaths_wide$normal)),]
fluxpaths_unique_zerofree <- unique(analyzed_fluxpaths_wide_noNANs$fluxpath_subsystem) #count the number of metabolites
#for this dataset, there are no nans and all data are kept
#replace remaining NA values with zeroes
analyzed_fluxpaths_wide_noNANs$tumor[is.na(analyzed_fluxpaths_wide_noNANs$tumor)] = 0
analyzed_fluxpaths_wide_noNANs$normal[is.na(analyzed_fluxpaths_wide_noNANs$normal)] = 0
#create a difference column
analyzed_fluxpaths_wide_noNANs$difference <- (analyzed_fluxpaths_wide_noNANs$tumor - analyzed_fluxpaths_wide_noNANs$normal)
#remove any rows where differences are 0
analyzed_fluxpaths_wide_noNANs_zerofree <- analyzed_fluxpaths_wide_noNANs[(analyzed_fluxpaths_wide_noNANs$difference !=0),]
#NOTE: for this dataset, a few rows are removed but no fluxpaths
#how many unique fluxes are we including here? 
fluxpaths_zerofree_unique <- unique(analyzed_fluxpaths_wide_noNANs_zerofree$fluxpath_subsystem)
#clean up
rm(fluxpaths_unique_zerofree, fluxpaths_zerofree_unique, analyzed_fluxpaths_wide_noNANs, analyzed_fluxpaths_wide, analyzed_fluxpaths_subset_averaged)

#########Step 3: filter by flux subsystems which are active in at least half the total number of samples 
#first, create a filtering dataframe that contains counts of each fluxpath_name incidence in the dataframe
#need to have n=2 if values are nonzero in both tumor and normal columns; n=1 if values are nonzero in only one column
analyzed_fluxpaths_filtering_df <- analyzed_fluxpaths_wide_noNANs_zerofree %>%  #create a new dataframe by filtering the old one
  mutate(count_of_nonzeros = rowSums(.[3:4]!=0)) %>% #create a new column, count_of_nonzeroes, that counts # nonzero values in columns 3 and 4 for each row
  group_by(fluxpath_subsystem) %>% #group by the column of interest
  summarize(n = sum(count_of_nonzeros)) %>% #sums the count_of_nonzeros values within each fluxpath_subsystem and saves the value as n
  mutate(n) #add n as a new column in the dataframe
#create a new dataframe that merges the old and new dataframes
analyzed_fluxpaths_filter_added <- dplyr::inner_join(analyzed_fluxpaths_wide_noNANs_zerofree, analyzed_fluxpaths_filtering_df, by= "fluxpath_subsystem")
#now, remove all rows where n (the name of the count column) is less than ...36, which is half the samples
analyzed_fluxpaths_filtered <- filter(analyzed_fluxpaths_filter_added, n >= 36)
#count number of unique metabolites left
fluxpaths_filtered_unique <- unique(analyzed_fluxpaths_filtered$fluxpath_subsystem)
#clean up
rm(fluxpaths_filtered_unique, analyzed_fluxpaths_wide_noNANs_zerofree)

##########Step 4: split into multiple dataframes by metabolite and do stats on each metabolite
#first, edit the dataframe
#drop a couple of columns
analyzed_fluxpaths_paired_only <- select(analyzed_fluxpaths_filtered, -c(difference, n))
#rename a couple of columns to make life easier later. Yeah, yeah, I know, this is the reverse of what I did before.
#analyzed_metabolites_paired_only <- analyzed_metabolites_paired_only %>% rename(normal= mean_genus_GR_normal, tumor= mean_genus_GR_tumor) #rename these columns
#convert to long form
analyzed_fluxpaths_paired_long <- reshape2::melt(data= analyzed_fluxpaths_paired_only,
                                                 id.vars= c("Patient_Blind_ID", "fluxpath_subsystem"),
                                                 variable.name = "Description",
                                                 value.name = "fluxpath_activity")
#split up by different subsystems and see where we get. 
split_analyzed_fluxpaths_by_flux <- split(analyzed_fluxpaths_paired_long, with(analyzed_fluxpaths_paired_long, interaction(fluxpath_subsystem)), drop = TRUE)
by_flux_stats <- lapply(split_analyzed_fluxpaths_by_flux, function(df){wilcox.test(fluxpath_activity~Description, data=df, exact= FALSE, paired= TRUE)})

##########Step 5: statistics
#make lists of p values
by_flux_pvals <- c() #initialize list
for (elem in by_flux_stats){
  new_value = elem$p.value
  by_flux_pvals <- c(by_flux_pvals, new_value)}
rm(new_value)
#make a list of q-values
by_flux_qvals <- p.adjust(by_flux_pvals, method = "fdr")
#make a list of genus_ids
flux_ids <- c()
for(elem in split_analyzed_fluxpaths_by_flux){
  new_value = elem$fluxpath_subsystem[1]
  flux_ids <- c(flux_ids, new_value)}
#merge all lists together. 
fluxpath_statistics <- data.frame(flux_ids, by_flux_pvals, by_flux_qvals)

###########Step 6: clean up and remove extra variables
rm(elem, split_analyzed_fluxpaths_by_flux, by_flux_pvals, by_flux_qvals, flux_ids, new_value)

###########Step 7: pull out significant p-values only
fluxpath_statistics <- fluxpath_statistics[order(fluxpath_statistics$by_flux_pvals),] #re-order by p-value

##########Step 8: save results as a csv file
fluxpath_statistics <- fluxpath_statistics[order(fluxpath_statistics$by_flux_qvals),] #re-order by q-value
fluxpath_statistics <- fluxpath_statistics %>% rename("fluxpath_name"="flux_ids")
write.csv(as.data.frame(fluxpath_statistics), file="01.10.22_Fusobacterium_flux_subsystems_stats_Niccolai2020_data.csv")

#######################################################ROUND 2:::: do all fluxes#############################################################

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Niccolai2020/") 

##########Step 2: import data and adjust factors
analyzed_fluxpaths <- read_csv("01.10.22_Niccolai2020_Fuso_fluxes_with_metadata_longform.csv") 
analyzed_fluxpaths$Patient_Blind_ID <- as.factor(analyzed_fluxpaths$Patient_Blind_ID)
analyzed_fluxpaths <- analyzed_fluxpaths %>% rename(fluxpath_name= fluxpath_id, 
                                                    fluxpath_subsystem= subsystem,
                                                    flux_value = fluxpath_values,
                                                    fluxpath_description = description)

#########Step 3a: adding across OTU_ids
#because we have multiple fuso_ids in this dataset, we need to collapse these. I don't want to fuck around with multiple OTUs per sample
#here, we would like to sum all flux_values across OTU_ids, within fluxpath_name/sample_ID combos
#first, drop unnecessary columns (i.e. what we're adding across)
sum_fuso_fluxpaths <- subset(analyzed_fluxpaths, select = -c(OTU_ID))
#create a new sorting column to get unique values for each Sample_ID/Fluxpath_name combo
sum_fuso_fluxpaths <- sum_fuso_fluxpaths %>% 
  mutate(sorting_col = paste0(SampleID, "_", fluxpath_name)) 
#Drop NAN values; converting to zeroes screws with averages
#sum_fuso_fluxpaths_noNANs <- sum_fuso_fluxpaths[!is.na(sum_fuso_fluxpaths$flux_value),]
#okay, here is another source of contention. I don't want to lose fluxpaths. And really, isn't zero basically NAN, since we've already
#filtered by Fuso-specific samples?
#NEW 12.21.21: replace NANs wtih zeroes
sum_fuso_fluxpaths$flux_value <- as.numeric(sum_fuso_fluxpaths$flux_value)
sum_fuso_fluxpaths$flux_value[is.na(sum_fuso_fluxpaths$flux_value)] <- 0

#then sum across sorting col:
sum_fuso_fluxpaths_noNANs_collapsed <- sum_fuso_fluxpaths %>%        # Specify data frame
  group_by(sorting_col, fluxpath_name, SampleID, fluxpath_description, Patient_Blind_ID, Description) %>%    # Specify group indicator (columns to keep)
  summarise_at(vars(flux_value),                   # Specify column
               list(sum_fluxpath_amount = sum)) %>% # Specify function
  ungroup() #do this so select (dropping columns) can happen later
#nice. End up with a dataframe where each row is a unique subsystem/sample average of fluxes from pathways within that subsystem
#how many unique fluxes are we including here? 
sum_fluxpaths_unique <- unique(sum_fuso_fluxpaths_noNANs_collapsed$fluxpath_name)

#rename for latter consistency in code
sum_fuso_fluxpaths_noNANs_collapsed <- sum_fuso_fluxpaths_noNANs_collapsed %>% rename(flux_value = sum_fluxpath_amount)

#clean up
rm(sum_fuso_fluxpaths, sum_fluxpaths_unique)

#########Step 4: collapsing by trace minerals: borrowed from 09.20.21_Burns2015_inputoutput_analysis_plus_traceminerals.R
#there are several metabolites (EX_ca2(e), EX_cobalt2(e), EX_cu2(e), EX_mg2(e), EX_mn2(e), EX_zn2(e), EX_k(e), EX_cl(e), EX_so4(e))
#which have identical fluxes in all communities. The goal is to collapse these into a single group, called "Trace minerals"
#first, create a vector of these metabolites
trace_minerals <- c("EX_ca2(e)", "EX_cobalt2(e)", "EX_cu2(e)", "EX_mg2(e)", "EX_mn2(e)", "EX_zn2(e)", "EX_k(e)", "EX_cl(e)", "EX_so4(e)")
#Figure out some way to count the number of trace minerals with nonzero values within each sample. 
#As long as it's always 10 or 0, can collapse
fuso_fluxpaths_only_filtered <- sum_fuso_fluxpaths_noNANs_collapsed

#first, create a dataframe that ONLY contains trace minerals
fuso_fluxpaths_trace_minerals_only <- fuso_fluxpaths_only_filtered[fuso_fluxpaths_only_filtered$fluxpath_name %in% trace_minerals, ]
#use dplyr to round to 6 significant digits.
fuso_fluxpaths_trace_minerals_only <- fuso_fluxpaths_trace_minerals_only %>% 
  mutate(flux_value=signif(flux_value, 6)) 
#then, create a filtering_DF that counts the number of unique values within each sample_ID_ish
trace_minerals_filtering_df <- fuso_fluxpaths_trace_minerals_only %>%  #create a new dataframe by filtering the old one
  group_by(SampleID) %>% #group by the column of interest
  summarize(n = n_distinct(flux_value)) #count the number of occurrences of each unique value in metabolite_amount and store as variable 'n'

#okay, now we have all 1s,This is actually good. It means within each sample_ID_ish,

#time to collapse! for each sample_ID_ish, only keep ONE of the rows where metabolite_names values are in the trace_minerals vector
#and rename metabolite_name to trace_mineral.
fuso_fluxpaths_collapsed <-data.frame(fuso_fluxpaths_only_filtered) #copy the original dataframe
#first, rename "EX_cobalt2(e)" to "exchange of trace elements"
fuso_fluxpaths_collapsed$fluxpath_name[fuso_fluxpaths_collapsed$fluxpath_name== "EX_cobalt2(e)"] <- "Trace_elements_exchange"
fuso_fluxpaths_collapsed$fluxpath_description[fuso_fluxpaths_collapsed$fluxpath_description== "Co2+ exchange"] <- "Exchange of nine trace elements"
#then, remove all rows where trace_minerals vector values are in the metabolite_name column
fuso_fluxpaths_collapsed<- fuso_fluxpaths_collapsed[!(fuso_fluxpaths_collapsed$fluxpath_name %in% trace_minerals),]
fuso_fluxpaths_collapsed_unique <- unique(fuso_fluxpaths_collapsed$fluxpath_name)
#clean up- note that we NEED fuso_fluxpaths_only_filtered
rm(fuso_fluxpaths_collapsed_unique, fuso_fluxpaths_trace_minerals_only, trace_minerals_filtering_df, trace_minerals, fuso_fluxpaths_noNANs_unique)

#########Step 5: remove samples where tumor and normal fluxes = 0
fuso_fluxpaths_wide <- select(fuso_fluxpaths_collapsed, -c(SampleID)) #drop SampleID column
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
#how many unique fluxes are we including here? 112
fuso_fluxpaths_zerofree_unique <- unique(fuso_fluxpaths_wide_noNANs_zerofree$fluxpath_name)
#cleanup: note to KEEP fuso_fluxpaths_collapsed
rm(fuso_fluxpaths_wide, fuso_fluxpaths_wide_noNANs, fuso_fluxpaths_noNANs_unique, fuso_fluxpaths_zerofree_unique)

########Step 6: filter by flux pathways that are active in ALL Fuso-containing samples
#first, find out how many samples have Fuso.
#EDIT here- should be using fuso_fluxpaths_only_filtered, not fuso_fluxpaths_collapsed (as the latter contains PBIDs
#which did NOT have fuso in both cases)
num_fuso_pos_samples <- unique(fuso_fluxpaths_only_filtered$SampleID) #72 samples, here.CHECK.

#first, create a filtering dataframe that contains counts of each fluxpath_name incidence in the dataframe
fuso_fluxpaths_filtering_df <- fuso_fluxpaths_wide_noNANs_zerofree %>%  #create a new dataframe by filtering the old one
  mutate(count_of_nonzeros = rowSums(.[5:6]!=0)) %>% #create a new column, count_of_nonzeroes, that counts # nonzero values in tumor and normal columns
  group_by(fluxpath_name) %>% #group by the column of interest
  summarize(n = sum(count_of_nonzeros)) %>% #sums the count_of_nonzeros values within each fluxpath_subsystem and saves the value as n
  mutate(n) #add n as a new column in the dataframe
#create a new dataframe that merges the old and new dataframes
fuso_fluxpaths_filter_added <- dplyr::inner_join(fuso_fluxpaths_wide_noNANs_zerofree, fuso_fluxpaths_filtering_df, by= "fluxpath_name")

#now, remove all rows where n (the name of the count column) is less than ... 44 is 100%, 22 is 50%, 12 is 20%
fuso_fluxpaths_abundance_filtered <- filter(fuso_fluxpaths_filter_added, n >= 36)
#count number of unique metabolites left: 43
fuso_fluxpaths_filtered_unique <- unique(fuso_fluxpaths_abundance_filtered$fluxpath_name)
#cleanup
rm(fuso_fluxpaths_filtering_df,fuso_fluxpaths_filter_added, fuso_fluxpaths_only_filtered, num_fuso_pos_samples, fuso_fluxpaths_filtered_unique, fuso_fluxpaths_wide_noNANs_zerofree)

##########Step 7: split into multiple dataframes by metabolite and do stats on each metabolite
#first, edit the dataframe; drop a couple of columns
#fuso_fluxpaths_abundance_filtered <- fuso_fluxpaths_abundance_filtered %>%  rename(fluxpath_description= description)
fuso_fluxpaths_paired_only <- select(fuso_fluxpaths_abundance_filtered, -c(difference, n, fluxpath_description, sorting_col))
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
fuso_fluxpath_statistics <- fuso_fluxpath_statistics %>% rename("fluxpath_name"="flux_ids")
#I would also like to add back the descriptions, please.
#fuso_fluxpaths_collapsed <- fuso_fluxpaths_collapsed %>%  rename(fluxpath_description= description)
fluxpath_key <- subset(fuso_fluxpaths_collapsed, select = c(fluxpath_name, fluxpath_description) )
fluxpath_key <- fluxpath_key %>% distinct(fluxpath_name, .keep_all = TRUE)
#create a new dataframe that merges the old and new dataframes
fuso_fluxpath_statistics2 <- dplyr::inner_join(fuso_fluxpath_statistics, fluxpath_key, by= "fluxpath_name")
write.csv(as.data.frame(fuso_fluxpath_statistics2), file="01.10.22_Fusobacterium_all_flux_stats_Niccolai2020_data.csv")

