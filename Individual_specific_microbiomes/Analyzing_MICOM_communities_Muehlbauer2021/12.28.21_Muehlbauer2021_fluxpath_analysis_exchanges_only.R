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

#starting point is 09.30.21_Hale2018_fluxpath_analysis_at_median_subsystem_level.R and 10.21.21_BUrns_fluxpath_analysis_exchanges_only.R
#and 10.26.21_Hale2018_fluxpath_analysis_exchanges_only.R

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Muehlbauer2021/") 

##########Step 1: import data and adjust factors
analyzed_fluxpaths <- read_csv("11.09.21_Muehlbauer2021_fluxes_with_metadata_longform.csv")
analyzed_fluxpaths$Library_name <- as.factor(analyzed_fluxpaths$Isolate)
analyzed_fluxpaths$experimental_treatment <- as.factor(analyzed_fluxpaths$experimental_treatment)
#colnames(analyzed_metabolites) #print column names as needed
analyzed_fluxpaths <- subset(analyzed_fluxpaths, select = -c(...1, fluxpath_formula, fluxpath_subsystem, Library_name))

#how many unique fluxes are we including here? 160
fluxpaths_unique <- unique(analyzed_fluxpaths$fluxpath_name) 

#########NEW Step 2c: collapse this dataset by subsystem across OTU_ids; here, SUM across OTUs within a fluxpath within a sample.
#drop unnecessary columns (i.e. what we're averating across)
sum_analyzed_fluxpaths_subset <- subset(analyzed_fluxpaths, select = -c(OTU_ID))
#create a new sorting column to get unique values for each Sample_ID/Fluxpath_subsystem combo
sum_analyzed_fluxpaths_subset <- sum_analyzed_fluxpaths_subset %>% 
  mutate(sorting_col = paste0(Sample_ID, "_", fluxpath_name)) 
#Drop NAN values; converting to zeroes screws with averages
sum_analyzed_fluxpaths_subset_noNANs <- sum_analyzed_fluxpaths_subset[!is.na(sum_analyzed_fluxpaths_subset$flux_value),]
#then sum across sorting col:
sum_analyzed_fluxpaths_subset_collapsed <- sum_analyzed_fluxpaths_subset_noNANs %>%        # Specify data frame
  group_by(sorting_col, Sample_ID, experimental_treatment, Isolate, fluxpath_name) %>%    # Specify group indicator (columns to keep)
  summarise_at(vars(flux_value),                   # Specify column
               list(sum_fluxpath_amount = sum)) %>% # Specify function
  ungroup() #do this so select (dropping columns) can happen later
#nice. End up with a dataframe where each row is a unique subsystem/sample average of fluxes from pathways within that subsystem
#how many unique fluxes are we including here? Still 160
sum_fluxpaths_unique <- unique(sum_analyzed_fluxpaths_subset_collapsed$fluxpath_name)

#######clean up
rm(sum_analyzed_fluxpaths_subset, sum_analyzed_fluxpaths_subset_noNANs)

#########Step 3c: filter to only include paired taxa (i.e. taxa which appear in both tumor and normal samples)
sum_analyzed_fluxpaths_wide <- select(sum_analyzed_fluxpaths_subset_collapsed, -c(sorting_col, Sample_ID)) #drop unnecessary columns
#first, convert to wide form
sum_analyzed_fluxpaths_wide <- sum_analyzed_fluxpaths_wide %>% 
  spread(experimental_treatment, sum_fluxpath_amount)
#then, drop all rows where both tumor and normal values are nan.
sum_analyzed_fluxpaths_wide_noNANs <- sum_analyzed_fluxpaths_wide[!(is.na(sum_analyzed_fluxpaths_wide$Colonocytes) & is.na(sum_analyzed_fluxpaths_wide$Only)
                                                                    & is.na(sum_analyzed_fluxpaths_wide$Uncultured)),]
sum_fluxpaths_unique_zerofree <- unique(sum_analyzed_fluxpaths_wide_noNANs$fluxpath_name) #count the number of metabolites
#replace remaining NA values with zeroes
sum_analyzed_fluxpaths_wide_noNANs$Colonocytes[is.na(sum_analyzed_fluxpaths_wide_noNANs$Colonocytes)] = 0
sum_analyzed_fluxpaths_wide_noNANs$Only[is.na(sum_analyzed_fluxpaths_wide_noNANs$Only)] = 0
sum_analyzed_fluxpaths_wide_noNANs$Uncultured[is.na(sum_analyzed_fluxpaths_wide_noNANs$Uncultured)] = 0
#create a difference column: NOT DOING IT HERE because there are three groups
#sum_analyzed_fluxpaths_wide_noNANs$difference <- (sum_analyzed_fluxpaths_wide_noNANs$tumor - sum_analyzed_fluxpaths_wide_noNANs$normal)
#remove any rows where differences are 0
#sum_analyzed_fluxpaths_wide_noNANs_zerofree <- sum_analyzed_fluxpaths_wide_noNANs[(sum_analyzed_fluxpaths_wide_noNANs$difference !=0),]
#NOTE: for this dataset, a few rows are removed but no fluxpaths
#how many unique fluxes are we including here? 
sum_fluxpaths_zerofree_unique <- unique(sum_analyzed_fluxpaths_wide_noNANs$fluxpath_name)
#cleanup
rm(sum_analyzed_fluxpaths_wide, sum_analyzed_fluxpaths_wide_noNANs, sum_fluxpaths_unique_zerofree)

#########Step 4c: filter by flux pathways that are active in ALL 88 samples
#first, create a filtering dataframe that contains counts of each fluxpath_name incidence in the dataframe
#need to have n=2 if values are nonzero in both tumor and normal columns; n=1 if values are nonzero in only one column
sum_analyzed_fluxpaths_filtering_df <- sum_analyzed_fluxpaths_wide_noNANs %>%  #create a new dataframe by filtering the old one
  mutate(count_of_nonzeros = rowSums(.[3:5]!=0)) %>% #create a new column, count_of_nonzeroes, that counts # nonzero values in tumor and normal columns
  group_by(fluxpath_name) %>% #group by the column of interest
  summarize(n = sum(count_of_nonzeros)) %>% #sums the count_of_nonzeros values within each fluxpath_subsystem and saves the value as n
  mutate(n) #add n as a new column in the dataframe
#create a new dataframe that merges the old and new dataframes
sum_analyzed_fluxpaths_filter_added <- dplyr::inner_join(sum_analyzed_fluxpaths_wide_noNANs, sum_analyzed_fluxpaths_filtering_df, by= "fluxpath_name")
#now, remove all rows where n (the name of the count column) is less than ...
sum_analyzed_fluxpaths_filtered <- filter(sum_analyzed_fluxpaths_filter_added, n >= 12)
#count number of unique metabolites left
sum_fluxpaths_filtered_unique <- unique(sum_analyzed_fluxpaths_filtered$fluxpath_name) #now we have 127
#cleanup
rm(sum_analyzed_fluxpaths_filtering_df)

##########Step 5c: split into multiple dataframes by metabolite and do stats on each metabolite
#first, edit the dataframe; drop a couple of columns
sum_analyzed_fluxpaths_paired_only <- select(sum_analyzed_fluxpaths_filtered, -c(n))
#convert to long form
sum_analyzed_fluxpaths_paired_long <- reshape2::melt(data= sum_analyzed_fluxpaths_paired_only,
                                                     id.vars= c("Isolate", "fluxpath_name"),
                                                     variable.name = "experimental_treatment",
                                                     value.name = "fluxpath_activity")

#SAVE HERE and move on. 
write.csv(as.data.frame(sum_analyzed_fluxpaths_paired_only), file="12.28.21_filtered_sum_exchange_fluxes_only_alldata_Muehlbauer2020_data.csv")


#split up by different subsystems and see where we get. 
sum_split_analyzed_fluxpaths_by_flux <- split(sum_analyzed_fluxpaths_paired_long, 
                                              with(sum_analyzed_fluxpaths_paired_long, interaction(fluxpath_name)), drop = TRUE)
sum_by_flux_stats <- lapply(sum_split_analyzed_fluxpaths_by_flux, function(df){wilcox.test(fluxpath_activity~experimental_treatment, data=df, exact= FALSE, paired= TRUE)})
#make lists of p values
by_flux_pvals <- c() #initialize list
for (elem in sum_by_flux_stats){
  new_value = elem$p.value
  by_flux_pvals <- c(by_flux_pvals, new_value)}
rm(new_value)
#make a list of q-values
by_flux_qvals <- p.adjust(by_flux_pvals, method = "fdr")
#make a list of genus_ids
flux_ids <- c()
for(elem in sum_split_analyzed_fluxpaths_by_flux){
  new_value = elem$fluxpath_name[1]
  flux_ids <- c(flux_ids, new_value)}
#merge all lists together. 
sum_fluxpath_statistics <- data.frame(flux_ids, by_flux_pvals, by_flux_qvals)
#clean up and remove extra variables
rm(elem, sum_split_analyzed_fluxpaths_by_flux, by_flux_pvals, by_flux_qvals, flux_ids, new_value, sum_by_flux_stats)

###########Step 6: sort and save
sum_fluxpath_statistics <- sum_fluxpath_statistics[order(sum_fluxpath_statistics$by_flux_qvals),] #re-order by q-value
sum_fluxpath_statistics <- sum_fluxpath_statistics %>% rename("fluxpath_name"="flux_ids" )
#I would also like to add back the Experimental_treatments, please.
fluxpath_key <- subset(analyzed_fluxpaths, select = c(fluxpath_name, fluxpath_Experimental_treatment) )
fluxpath_key <- fluxpath_key %>% distinct(fluxpath_name, .keep_all = TRUE)
#create a new dataframe that merges the old and new dataframes
sum_fluxpath_statistics2 <- dplyr::inner_join(sum_fluxpath_statistics, fluxpath_key, by= "fluxpath_name")
write.csv(as.data.frame(sum_fluxpath_statistics2), file="11.09.21_filtered_sum_exchange_fluxes_only_stats_Niccolai2020_data.csv")

##########Step 7: prepare data for graphing 
#add a new column called diff. tumor-normal
sum_analyzed_fluxpaths_paired_only$diff <- (sum_analyzed_fluxpaths_paired_only$tumor - sum_analyzed_fluxpaths_paired_only$normal)
#NOW we can add a new column with useful binning of differences for graphing purposes
sum_analyzed_fluxpaths_paired_wide <- sum_analyzed_fluxpaths_paired_only %>%
  mutate(Difference_direction = case_when(
    diff > 0 ~ "tumor > normal",
    diff < 0 ~ "tumor < normal",
    diff == 0 ~ "Zero change",))
#I would also like to add back the Experimental_treatments, please. (The below two lines should already be done, from above)
#fluxpath_key <- subset(fuso_fluxpaths_collapsed, select = c(fluxpath_name, fluxpath_Experimental_treatment) )
#fluxpath_key <- fluxpath_key %>% distinct(fluxpath_name, .keep_all = TRUE)
#create a new dataframe that merges the old and new dataframes
sum_analyzed_fluxpaths_paried_wide2 <- dplyr::inner_join(sum_analyzed_fluxpaths_paired_wide, fluxpath_key, by= "fluxpath_name")
write.csv(as.data.frame(sum_analyzed_fluxpaths_paried_wide2), file="11.09.21_filtered_sum_exchange_fluxes_only_alldata_Muehlbauer2020_data.csv")

#NOW, graph.

#Exchange of deoxycytidine
ggpaired(subset(sum_analyzed_fluxpaths_paired_wide, fluxpath_name %in% c("EX_dcyt(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Exchange of deoxycytidine\n sum of fluxes, Niccolai2020 data"),
       subtitle = expression("paired Wilcoxan, p=0.00178, q= 0.0630")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Sum of fluxes,\nmmol/(gDW * h)", limits=c(-10000,5000)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Exchange of riboflavin
ggpaired(subset(sum_analyzed_fluxpaths_paired_wide, fluxpath_name %in% c("EX_ribflv(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Exchange of riboflavin\n sum of fluxes, Niccolai2020 data"),
       subtitle = expression("paired Wilcoxan, p=0.00112, q= 0.0630")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Sum of fluxes,\nmmol/(gDW * h)", limits=c(0,8000)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Exchange of Xanthine
ggpaired(subset(sum_analyzed_fluxpaths_paired_wide, fluxpath_name %in% c("EX_xan(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Exchange of xanthine\n sum of fluxes, Niccolai2020 data"),
       subtitle = expression("paired Wilcoxan, p=0.00186, q= 0.0630")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Sum of fluxes,\nmmol/(gDW * h)", limits=c(0,30000)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Exchange of Uracil
ggpaired(subset(sum_analyzed_fluxpaths_paired_wide, fluxpath_name %in% c("EX_ura(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Exchange of uracil\n sum of fluxes, Niccolai2020 data"),
       subtitle = expression("paired Wilcoxan, p=0.00223, q= 0.0630")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Sum of fluxes,\nmmol/(gDW * h)", limits=c(0,50000)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")
