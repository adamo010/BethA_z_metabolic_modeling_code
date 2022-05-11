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

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Hale2018/") 

##########Step 1: import data and adjust factors
analyzed_fluxpaths <- read_csv("10.13.21_Hale2018_fluxes_with_metadata_longform.csv")
analyzed_fluxpaths$host_subject_id <- as.factor(analyzed_fluxpaths$host_subject_id)
analyzed_fluxpaths$normal_adjacent_or_tumor_tissue_specimen <- as.factor(analyzed_fluxpaths$normal_adjacent_or_tumor_tissue_specimen)
#colnames(analyzed_metabolites) #print column names as needed
analyzed_fluxpaths <- subset(analyzed_fluxpaths, select = -c(...1, fluxpath_formula, fluxpath_subsystem, sampleIDs_from_read_counts))

#how many unique fluxes are we including here? 
fluxpaths_unique <- unique(analyzed_fluxpaths$fluxpath_name)

#edit column names to match Burns because I am lazy and don't want to think too much about rerunning this code again.
analyzed_fluxpaths <- analyzed_fluxpaths %>% rename(Patient_Blind_ID = host_subject_id, Description = normal_adjacent_or_tumor_tissue_specimen)

#########NEW Step 2a: collapse this dataset by subsystem across OTU_ids; here, average across OTUs within a fluxpath within a sample.
#drop unnecessary columns (i.e. what we're averating across)
av_analyzed_fluxpaths_subset <- subset(analyzed_fluxpaths, select = -c(OTU_ID))
#create a new sorting column to get unique values for each Sample_ID/Fluxpath_subsystem combo
av_analyzed_fluxpaths_subset <- av_analyzed_fluxpaths_subset %>% 
  mutate(sorting_col = paste0(Sample_ID, "_", fluxpath_name)) 
#Drop NAN values; converting to zeroes screws with averages
av_analyzed_fluxpaths_subset_noNANs <- av_analyzed_fluxpaths_subset[!is.na(av_analyzed_fluxpaths_subset$flux_value),]
#NEW 09.29.21: take median instead of mean, to control for wild -ve values in some paths. 
#then average across sorting col:
av_analyzed_fluxpaths_subset_collapsed <- av_analyzed_fluxpaths_subset_noNANs %>%        # Specify data frame
  group_by(sorting_col, Sample_ID, Description, Patient_Blind_ID, fluxpath_name) %>%    # Specify group indicator (columns to keep)
  summarise_at(vars(flux_value),                   # Specify column
               list(av_fluxpath_amount = mean)) %>% # Specify function
  ungroup() #do this so select (dropping columns) can happen later
#nice. End up with a dataframe where each row is a unique subsystem/sample average of fluxes from pathways within that subsystem
#how many unique fluxes are we including here? 
av_fluxpaths_unique <- unique(av_analyzed_fluxpaths_subset_collapsed$fluxpath_name)

#########NEW Step 2b: collapse this dataset by subsystem across OTU_ids; here, median value across OTUs within a fluxpath within a sample.
#drop unnecessary columns (i.e. what we're averating across)
med_analyzed_fluxpaths_subset <- subset(analyzed_fluxpaths, select = -c(OTU_ID))
#create a new sorting column to get unique values for each Sample_ID/Fluxpath_subsystem combo
med_analyzed_fluxpaths_subset <- med_analyzed_fluxpaths_subset %>% 
  mutate(sorting_col = paste0(Sample_ID, "_", fluxpath_name)) 
#Drop NAN values; converting to zeroes screws with averages
med_analyzed_fluxpaths_subset_noNANs <- med_analyzed_fluxpaths_subset[!is.na(med_analyzed_fluxpaths_subset$flux_value),]
#NEW 09.29.21: take median instead of mean, to control for wild -ve values in some paths. 
#then average across sorting col:
med_analyzed_fluxpaths_subset_collapsed <- med_analyzed_fluxpaths_subset_noNANs %>%        # Specify data frame
  group_by(sorting_col, Sample_ID, Description, Patient_Blind_ID, fluxpath_name) %>%    # Specify group indicator (columns to keep)
  summarise_at(vars(flux_value),                   # Specify column
               list(med_fluxpath_amount = median)) %>% # Specify function
  ungroup() #do this so select (dropping columns) can happen later
#nice. End up with a dataframe where each row is a unique subsystem/sample average of fluxes from pathways within that subsystem
#how many unique fluxes are we including here? 
med_fluxpaths_unique <- unique(med_analyzed_fluxpaths_subset_collapsed$fluxpath_name)

#########NEW Step 2c: collapse this dataset by subsystem across OTU_ids; here, SUM across OTUs within a fluxpath within a sample.
#drop unnecessary columns (i.e. what we're averating across)
sum_analyzed_fluxpaths_subset <- subset(analyzed_fluxpaths, select = -c(OTU_ID))
#create a new sorting column to get unique values for each Sample_ID/Fluxpath_subsystem combo
sum_analyzed_fluxpaths_subset <- sum_analyzed_fluxpaths_subset %>% 
  mutate(sorting_col = paste0(Sample_ID, "_", fluxpath_name)) 
#Drop NAN values; converting to zeroes screws with averages
sum_analyzed_fluxpaths_subset_noNANs <- sum_analyzed_fluxpaths_subset[!is.na(sum_analyzed_fluxpaths_subset$flux_value),]
#NEW 09.29.21: take median instead of mean, to control for wild -ve values in some paths. 
#then average across sorting col:
sum_analyzed_fluxpaths_subset_collapsed <- sum_analyzed_fluxpaths_subset_noNANs %>%        # Specify data frame
  group_by(sorting_col, Sample_ID, Description, Patient_Blind_ID, fluxpath_name) %>%    # Specify group indicator (columns to keep)
  summarise_at(vars(flux_value),                   # Specify column
               list(sum_fluxpath_amount = sum)) %>% # Specify function
  ungroup() #do this so select (dropping columns) can happen later
#nice. End up with a dataframe where each row is a unique subsystem/sample average of fluxes from pathways within that subsystem
#how many unique fluxes are we including here? 
sum_fluxpaths_unique <- unique(sum_analyzed_fluxpaths_subset_collapsed$fluxpath_name)

#######clean up
rm(av_analyzed_fluxpaths_subset, av_analyzed_fluxpaths_subset_noNANs, med_analyzed_fluxpaths_subset, med_analyzed_fluxpaths_subset_noNANs,
   sum_analyzed_fluxpaths_subset, sum_analyzed_fluxpaths_subset_noNANs)

#########Step 3a: filter to only include paired taxa (i.e. taxa which appear in both tumor and normal samples)
av_analyzed_fluxpaths_wide <- select(av_analyzed_fluxpaths_subset_collapsed, -c(sorting_col, Sample_ID)) #drop unnecessary columns
#first, convert to wide form
av_analyzed_fluxpaths_wide <- av_analyzed_fluxpaths_wide %>% 
  spread(Description, av_fluxpath_amount)
#then, drop all rows where both tumor and normal values are nan.
av_analyzed_fluxpaths_wide_noNANs <- av_analyzed_fluxpaths_wide[!(is.na(av_analyzed_fluxpaths_wide$tumor) & is.na(av_analyzed_fluxpaths_wide$normal)),]
av_fluxpaths_unique_zerofree <- unique(av_analyzed_fluxpaths_wide_noNANs$fluxpath_name) #count the number of metabolites
#replace remaining NA values with zeroes
av_analyzed_fluxpaths_wide_noNANs$tumor[is.na(av_analyzed_fluxpaths_wide_noNANs$tumor)] = 0
av_analyzed_fluxpaths_wide_noNANs$normal[is.na(av_analyzed_fluxpaths_wide_noNANs$normal)] = 0
#create a difference column
av_analyzed_fluxpaths_wide_noNANs$difference <- (av_analyzed_fluxpaths_wide_noNANs$tumor - av_analyzed_fluxpaths_wide_noNANs$normal)
#remove any rows where differences are 0
av_analyzed_fluxpaths_wide_noNANs_zerofree <- av_analyzed_fluxpaths_wide_noNANs[(av_analyzed_fluxpaths_wide_noNANs$difference !=0),]
#NOTE: for this dataset, a few rows are removed but no fluxpaths
#how many unique fluxes are we including here? 
av_fluxpaths_zerofree_unique <- unique(av_analyzed_fluxpaths_wide_noNANs_zerofree$fluxpath_name)
#cleanup
rm(av_analyzed_fluxpaths_wide, av_analyzed_fluxpaths_wide_noNANs, av_fluxpaths_unique_zerofree)

#########Step 3b: filter to only include paired taxa (i.e. taxa which appear in both tumor and normal samples)
med_analyzed_fluxpaths_wide <- select(med_analyzed_fluxpaths_subset_collapsed, -c(sorting_col, Sample_ID)) #drop unnecessary columns
#first, convert to wide form
med_analyzed_fluxpaths_wide <- med_analyzed_fluxpaths_wide %>% 
  spread(Description, med_fluxpath_amount)
#then, drop all rows where both tumor and normal values are nan.
med_analyzed_fluxpaths_wide_noNANs <- med_analyzed_fluxpaths_wide[!(is.na(med_analyzed_fluxpaths_wide$tumor) & is.na(med_analyzed_fluxpaths_wide$normal)),]
med_fluxpaths_unique_zerofree <- unique(med_analyzed_fluxpaths_wide_noNANs$fluxpath_name) #count the number of metabolites
#replace remaining NA values with zeroes
med_analyzed_fluxpaths_wide_noNANs$tumor[is.na(med_analyzed_fluxpaths_wide_noNANs$tumor)] = 0
med_analyzed_fluxpaths_wide_noNANs$normal[is.na(med_analyzed_fluxpaths_wide_noNANs$normal)] = 0
#create a difference column
med_analyzed_fluxpaths_wide_noNANs$difference <- (med_analyzed_fluxpaths_wide_noNANs$tumor - med_analyzed_fluxpaths_wide_noNANs$normal)
#remove any rows where differences are 0
med_analyzed_fluxpaths_wide_noNANs_zerofree <- med_analyzed_fluxpaths_wide_noNANs[(med_analyzed_fluxpaths_wide_noNANs$difference !=0),]
#NOTE: for this dataset, a few rows are removed but no fluxpaths
#how many unique fluxes are we including here? 
med_fluxpaths_zerofree_unique <- unique(med_analyzed_fluxpaths_wide_noNANs_zerofree$fluxpath_name)
#cleanup
rm(med_analyzed_fluxpaths_wide, med_analyzed_fluxpaths_wide_noNANs, med_fluxpaths_unique_zerofree)

#########Step 3c: filter to only include paired taxa (i.e. taxa which appear in both tumor and normal samples)
sum_analyzed_fluxpaths_wide <- select(sum_analyzed_fluxpaths_subset_collapsed, -c(sorting_col, Sample_ID)) #drop unnecessary columns
#first, convert to wide form
sum_analyzed_fluxpaths_wide <- sum_analyzed_fluxpaths_wide %>% 
  spread(Description, sum_fluxpath_amount)
#then, drop all rows where both tumor and normal values are nan.
sum_analyzed_fluxpaths_wide_noNANs <- sum_analyzed_fluxpaths_wide[!(is.na(sum_analyzed_fluxpaths_wide$tumor) & is.na(sum_analyzed_fluxpaths_wide$normal)),]
sum_fluxpaths_unique_zerofree <- unique(sum_analyzed_fluxpaths_wide_noNANs$fluxpath_name) #count the number of metabolites
#replace remaining NA values with zeroes
sum_analyzed_fluxpaths_wide_noNANs$tumor[is.na(sum_analyzed_fluxpaths_wide_noNANs$tumor)] = 0
sum_analyzed_fluxpaths_wide_noNANs$normal[is.na(sum_analyzed_fluxpaths_wide_noNANs$normal)] = 0
#create a difference column
sum_analyzed_fluxpaths_wide_noNANs$difference <- (sum_analyzed_fluxpaths_wide_noNANs$tumor - sum_analyzed_fluxpaths_wide_noNANs$normal)
#remove any rows where differences are 0
sum_analyzed_fluxpaths_wide_noNANs_zerofree <- sum_analyzed_fluxpaths_wide_noNANs[(sum_analyzed_fluxpaths_wide_noNANs$difference !=0),]
#NOTE: for this dataset, a few rows are removed but no fluxpaths
#how many unique fluxes are we including here? 
sum_fluxpaths_zerofree_unique <- unique(sum_analyzed_fluxpaths_wide_noNANs_zerofree$fluxpath_name)
#cleanup
rm(sum_analyzed_fluxpaths_wide, sum_analyzed_fluxpaths_wide_noNANs, sum_fluxpaths_unique_zerofree)

#########Step 4a: filter by flux pathways that are active in ALL 88 samples
#first, create a filtering dataframe that contains counts of each fluxpath_name incidence in the dataframe
#need to have n=2 if values are nonzero in both tumor and normal columns; n=1 if values are nonzero in only one column
av_analyzed_fluxpaths_filtering_df <- av_analyzed_fluxpaths_wide_noNANs_zerofree %>%  #create a new dataframe by filtering the old one
  mutate(count_of_nonzeros = rowSums(.[3:4]!=0)) %>% #create a new column, count_of_nonzeroes, that counts # nonzero values in columns 3 and 4 for each row
  group_by(fluxpath_name) %>% #group by the column of interest
  summarize(n = sum(count_of_nonzeros)) %>% #sums the count_of_nonzeros values within each fluxpath_subsystem and saves the value as n
  mutate(n) #add n as a new column in the dataframe
#create a new dataframe that merges the old and new dataframes
av_analyzed_fluxpaths_filter_added <- dplyr::inner_join(av_analyzed_fluxpaths_wide_noNANs_zerofree, av_analyzed_fluxpaths_filtering_df, by= "fluxpath_name")
#now, remove all rows where n (the name of the count column) is less than ...
av_analyzed_fluxpaths_filtered <- filter(av_analyzed_fluxpaths_filter_added, n >= 88)
#count number of unique metabolites left
av_fluxpaths_filtered_unique <- unique(av_analyzed_fluxpaths_filtered$fluxpath_name)
#cleanup
rm(av_analyzed_fluxpaths_filtering_df)

#########Step 4b: filter by flux pathways that are active in ALL 88 samples
#first, create a filtering dataframe that contains counts of each fluxpath_name incidence in the dataframe
#need to have n=2 if values are nonzero in both tumor and normal columns; n=1 if values are nonzero in only one column
med_analyzed_fluxpaths_filtering_df <- med_analyzed_fluxpaths_wide_noNANs_zerofree %>%  #create a new dataframe by filtering the old one
  mutate(count_of_nonzeros = rowSums(.[3:4]!=0)) %>% #create a new column, count_of_nonzeroes, that counts # nonzero values in columns 3 and 4 for each row
  group_by(fluxpath_name) %>% #group by the column of interest
  summarize(n = sum(count_of_nonzeros)) %>% #sums the count_of_nonzeros values within each fluxpath_subsystem and saves the value as n
  mutate(n) #add n as a new column in the dataframe
#create a new dataframe that merges the old and new dataframes
med_analyzed_fluxpaths_filter_added <- dplyr::inner_join(med_analyzed_fluxpaths_wide_noNANs_zerofree, med_analyzed_fluxpaths_filtering_df, by= "fluxpath_name")
#now, remove all rows where n (the name of the count column) is less than ...
med_analyzed_fluxpaths_filtered <- filter(med_analyzed_fluxpaths_filter_added, n >= 88)
#count number of unique metabolites left
med_fluxpaths_filtered_unique <- unique(med_analyzed_fluxpaths_filtered$fluxpath_name)
#cleanup
rm(med_analyzed_fluxpaths_filtering_df)

#########Step 4c: filter by flux pathways that are active in ALL 88 samples
#first, create a filtering dataframe that contains counts of each fluxpath_name incidence in the dataframe
#need to have n=2 if values are nonzero in both tumor and normal columns; n=1 if values are nonzero in only one column
sum_analyzed_fluxpaths_filtering_df <- sum_analyzed_fluxpaths_wide_noNANs_zerofree %>%  #create a new dataframe by filtering the old one
  mutate(count_of_nonzeros = rowSums(.[3:4]!=0)) %>% #create a new column, count_of_nonzeroes, that counts # nonzero values in columns 3 and 4 for each row
  group_by(fluxpath_name) %>% #group by the column of interest
  summarize(n = sum(count_of_nonzeros)) %>% #sums the count_of_nonzeros values within each fluxpath_subsystem and saves the value as n
  mutate(n) #add n as a new column in the dataframe
#create a new dataframe that merges the old and new dataframes
sum_analyzed_fluxpaths_filter_added <- dplyr::inner_join(sum_analyzed_fluxpaths_wide_noNANs_zerofree, sum_analyzed_fluxpaths_filtering_df, by= "fluxpath_name")
#now, remove all rows where n (the name of the count column) is less than ...
sum_analyzed_fluxpaths_filtered <- filter(sum_analyzed_fluxpaths_filter_added, n >= 88)
#count number of unique metabolites left
sum_fluxpaths_filtered_unique <- unique(sum_analyzed_fluxpaths_filtered$fluxpath_name)
#cleanup
rm(sum_analyzed_fluxpaths_filtering_df)

##########Step 5a: split into multiple dataframes by metabolite and do stats on each metabolite
#first, edit the dataframe; drop a couple of columns
av_analyzed_fluxpaths_paired_only <- select(av_analyzed_fluxpaths_filtered, -c(difference, n))
#convert to long form
av_analyzed_fluxpaths_paired_long <- reshape2::melt(data= av_analyzed_fluxpaths_paired_only,
                                                    id.vars= c("Patient_Blind_ID", "fluxpath_name"),
                                                    variable.name = "Description",
                                                    value.name = "fluxpath_activity")
#split up by different subsystems and see where we get. 
av_split_analyzed_fluxpaths_by_flux <- split(av_analyzed_fluxpaths_paired_long, 
                                             with(av_analyzed_fluxpaths_paired_long, interaction(fluxpath_name)), drop = TRUE)
av_by_flux_stats <- lapply(av_split_analyzed_fluxpaths_by_flux, function(df){wilcox.test(fluxpath_activity~Description, data=df, exact= FALSE, paired= TRUE)})
#make lists of p values
by_flux_pvals <- c() #initialize list
for (elem in av_by_flux_stats){
  new_value = elem$p.value
  by_flux_pvals <- c(by_flux_pvals, new_value)}
rm(new_value)
#make a list of q-values
by_flux_qvals <- p.adjust(by_flux_pvals, method = "fdr")
#make a list of genus_ids
flux_ids <- c()
for(elem in av_split_analyzed_fluxpaths_by_flux){
  new_value = elem$fluxpath_name[1]
  flux_ids <- c(flux_ids, new_value)}
#merge all lists together. 
av_fluxpath_statistics <- data.frame(flux_ids, by_flux_pvals, by_flux_qvals)
#clean up and remove extra variables
rm(elem, av_split_analyzed_fluxpaths_by_flux, by_flux_pvals, by_flux_qvals, flux_ids, new_value, av_by_flux_stats)

##########Step 5b: split into multiple dataframes by metabolite and do stats on each metabolite
#first, edit the dataframe; drop a couple of columns
med_analyzed_fluxpaths_paired_only <- select(med_analyzed_fluxpaths_filtered, -c(difference, n))
#convert to long form
med_analyzed_fluxpaths_paired_long <- reshape2::melt(data= med_analyzed_fluxpaths_paired_only,
                                                     id.vars= c("Patient_Blind_ID", "fluxpath_name"),
                                                     variable.name = "Description",
                                                     value.name = "fluxpath_activity")
#split up by different subsystems and see where we get. 
med_split_analyzed_fluxpaths_by_flux <- split(med_analyzed_fluxpaths_paired_long, 
                                              with(med_analyzed_fluxpaths_paired_long, interaction(fluxpath_name)), drop = TRUE)
med_by_flux_stats <- lapply(med_split_analyzed_fluxpaths_by_flux, function(df){wilcox.test(fluxpath_activity~Description, data=df, exact= FALSE, paired= TRUE)})
#make lists of p values
by_flux_pvals <- c() #initialize list
for (elem in med_by_flux_stats){
  new_value = elem$p.value
  by_flux_pvals <- c(by_flux_pvals, new_value)}
rm(new_value)
#make a list of q-values
by_flux_qvals <- p.adjust(by_flux_pvals, method = "fdr")
#make a list of genus_ids
flux_ids <- c()
for(elem in med_split_analyzed_fluxpaths_by_flux){
  new_value = elem$fluxpath_name[1]
  flux_ids <- c(flux_ids, new_value)}
#merge all lists together. 
med_fluxpath_statistics <- data.frame(flux_ids, by_flux_pvals, by_flux_qvals)
#clean up and remove extra variables
rm(elem, med_split_analyzed_fluxpaths_by_flux, by_flux_pvals, by_flux_qvals, flux_ids, new_value, med_by_flux_stats)

##########Step 5c: split into multiple dataframes by metabolite and do stats on each metabolite
#first, edit the dataframe; drop a couple of columns
sum_analyzed_fluxpaths_paired_only <- select(sum_analyzed_fluxpaths_filtered, -c(difference, n))
#convert to long form
sum_analyzed_fluxpaths_paired_long <- reshape2::melt(data= sum_analyzed_fluxpaths_paired_only,
                                                     id.vars= c("Patient_Blind_ID", "fluxpath_name"),
                                                     variable.name = "Description",
                                                     value.name = "fluxpath_activity")
#split up by different subsystems and see where we get. 
sum_split_analyzed_fluxpaths_by_flux <- split(sum_analyzed_fluxpaths_paired_long, 
                                              with(sum_analyzed_fluxpaths_paired_long, interaction(fluxpath_name)), drop = TRUE)
sum_by_flux_stats <- lapply(sum_split_analyzed_fluxpaths_by_flux, function(df){wilcox.test(fluxpath_activity~Description, data=df, exact= FALSE, paired= TRUE)})
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
av_fluxpath_statistics <- av_fluxpath_statistics[order(av_fluxpath_statistics$by_flux_pvals),] #re-order by p-value
med_fluxpath_statistics <- med_fluxpath_statistics[order(med_fluxpath_statistics$by_flux_pvals),] #re-order by p-value
sum_fluxpath_statistics <- sum_fluxpath_statistics[order(sum_fluxpath_statistics$by_flux_pvals),] #re-order by p-value

write.csv(as.data.frame(av_fluxpath_statistics), file="10.26.21_filtered_mean_exchange_fluxes_only_stats_Hale2018_data.csv")
write.csv(as.data.frame(med_fluxpath_statistics), file="10.26.21_filtered_median_exchange_fluxes_only_stats_Hale2018_data.csv")
write.csv(as.data.frame(sum_fluxpath_statistics), file="10.26.21_filtered_sum_exchange_fluxes_only_stats_Hale2018_data.csv")

##########Step 7: prepare data for graphing 
#add a new column called diff. tumor-normal
av_analyzed_fluxpaths_paired_only$diff <- (av_analyzed_fluxpaths_paired_only$tumor - av_analyzed_fluxpaths_paired_only$normal)
#NOW we can add a new column with useful binning of differences for graphing purposes
av_analyzed_fluxpaths_paired_wide <- av_analyzed_fluxpaths_paired_only %>%
  mutate(Difference_direction = case_when(
    diff > 0 ~ "tumor > normal",
    diff < 0 ~ "tumor < normal",
    diff == 0 ~ "Zero change",))
write.csv(as.data.frame(av_analyzed_fluxpaths_paired_wide), file="10.26.21_filtered_mean_exchange_fluxes_only_alldata_Hale2018_data.csv")
#add a new column called diff. tumor-normal
med_analyzed_fluxpaths_paired_only$diff <- (med_analyzed_fluxpaths_paired_only$tumor - med_analyzed_fluxpaths_paired_only$normal)
#NOW we can add a new column with useful binning of differences for graphing purposes
med_analyzed_fluxpaths_paired_wide <- med_analyzed_fluxpaths_paired_only %>%
  mutate(Difference_direction = case_when(
    diff > 0 ~ "tumor > normal",
    diff < 0 ~ "tumor < normal",
    diff == 0 ~ "Zero change",))
write.csv(as.data.frame(med_analyzed_fluxpaths_paired_wide), file="10.26.21_filtered_median_exchange_fluxes_only_alldata_Hale2018_data.csv")
#add a new column called diff. tumor-normal
sum_analyzed_fluxpaths_paired_only$diff <- (sum_analyzed_fluxpaths_paired_only$tumor - sum_analyzed_fluxpaths_paired_only$normal)
#NOW we can add a new column with useful binning of differences for graphing purposes
sum_analyzed_fluxpaths_paired_wide <- sum_analyzed_fluxpaths_paired_only %>%
  mutate(Difference_direction = case_when(
    diff > 0 ~ "tumor > normal",
    diff < 0 ~ "tumor < normal",
    diff == 0 ~ "Zero change",))
write.csv(as.data.frame(sum_analyzed_fluxpaths_paired_wide), file="10.26.21_filtered_sum_exchange_fluxes_only_alldata_Hale2018_data.csv")

#NOW, graph.
#what are we graphing? Exchange of L-Histidine (EX_his_L(e))for both mean and sum

#Exchange of L-Histidine, mean
ggpaired(subset(av_analyzed_fluxpaths_paired_wide, fluxpath_name %in% c("EX_his_L(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Exchange of L-Histidine\n mean flux, Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.0148, q= 0.720")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Mean flux,\nmmol/(gDW * h)", limits=c(-80,0)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Exchange of L-Histidine, sum
ggpaired(subset(sum_analyzed_fluxpaths_paired_wide, fluxpath_name %in% c("EX_his_L(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Exchange of L-Histidine\n sum of fluxes, Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.0281, q= 0.594")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Sum of fluxes,\nmmol/(gDW * h)", limits=c(-6000,0)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Exchange of ammonia, mean
ggpaired(subset(av_analyzed_fluxpaths_paired_wide, fluxpath_name %in% c("EX_nh4(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Exchange of Ammonia\n mean flux, Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.0148, q= 0.720")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Mean flux,\nmmol/(gDW * h)", limits=c(0,600)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Exchange of ammonia, median
ggpaired(subset(med_analyzed_fluxpaths_paired_wide, fluxpath_name %in% c("EX_nh4(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Exchange of Ammonia\n median flux, Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.0617, q= 0.735")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Median flux,\nmmol/(gDW * h)", limits=c(0,700)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")


