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

#starting point is 09.29.21_Burns2015_fluxpath_analysis_at_median_subsystem_level.R

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/") 

##########Step 1: import data and adjust factors
analyzed_fluxpaths <- read_csv("10.13.21_Burns2015_fluxes_with_metadata_longform.csv")
analyzed_fluxpaths$Patient_Blind_ID <- as.factor(analyzed_fluxpaths$Patient_Blind_ID)
#colnames(analyzed_metabolites) #print column names as needed
analyzed_fluxpaths <- subset(analyzed_fluxpaths, select = -c(...1, fluxpath_formula, fluxpath_subsystem))

#how many unique fluxes are we including here? 
fluxpaths_unique <- unique(analyzed_fluxpaths$fluxpath_name)

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

write.csv(as.data.frame(av_fluxpath_statistics), file="10.25.21_filtered_mean_exchange_fluxes_only_stats_Burns2015_data.csv")
write.csv(as.data.frame(med_fluxpath_statistics), file="10.25.21_filtered_median_exchange_fluxes_only_stats_Burns2015_data.csv")
write.csv(as.data.frame(sum_fluxpath_statistics), file="10.25.21_filtered_sum_exchange_fluxes_only_stats_Burns2015_data.csv")

##########Step 7: prepare data for graphing 
#add a new column called diff. tumor-normal
av_analyzed_fluxpaths_paired_only$diff <- (av_analyzed_fluxpaths_paired_only$tumor - av_analyzed_fluxpaths_paired_only$normal)
#NOW we can add a new column with useful binning of differences for graphing purposes
av_analyzed_fluxpaths_paired_wide <- av_analyzed_fluxpaths_paired_only %>%
  mutate(Difference_direction = case_when(
    diff > 0 ~ "tumor > normal",
    diff < 0 ~ "tumor < normal",
    diff == 0 ~ "Zero change",))
write.csv(as.data.frame(av_analyzed_fluxpaths_paired_wide), file="10.25.21_filtered_mean_exchange_fluxes_only_alldata_Burns2015_data.csv")
#add a new column called diff. tumor-normal
med_analyzed_fluxpaths_paired_only$diff <- (med_analyzed_fluxpaths_paired_only$tumor - med_analyzed_fluxpaths_paired_only$normal)
#NOW we can add a new column with useful binning of differences for graphing purposes
med_analyzed_fluxpaths_paired_wide <- med_analyzed_fluxpaths_paired_only %>%
  mutate(Difference_direction = case_when(
    diff > 0 ~ "tumor > normal",
    diff < 0 ~ "tumor < normal",
    diff == 0 ~ "Zero change",))
write.csv(as.data.frame(med_analyzed_fluxpaths_paired_wide), file="10.25.21_filtered_median_exchange_fluxes_only_alldata_Burns2015_data.csv")
#add a new column called diff. tumor-normal
sum_analyzed_fluxpaths_paired_only$diff <- (sum_analyzed_fluxpaths_paired_only$tumor - sum_analyzed_fluxpaths_paired_only$normal)
#NOW we can add a new column with useful binning of differences for graphing purposes
sum_analyzed_fluxpaths_paired_wide <- sum_analyzed_fluxpaths_paired_only %>%
  mutate(Difference_direction = case_when(
    diff > 0 ~ "tumor > normal",
    diff < 0 ~ "tumor < normal",
    diff == 0 ~ "Zero change",))
write.csv(as.data.frame(sum_analyzed_fluxpaths_paired_wide), file="10.25.21_filtered_sum_exchange_fluxes_only_alldata_Burns2015_data.csv")

#NOW, graph.
#what are we graphing? Exchange of L-Asparagine (EX_asn_L(e))for both mean and sum

#Exchange of L-Asparagine, mean
ggpaired(subset(av_analyzed_fluxpaths_paired_wide, fluxpath_name %in% c("EX_asn_L(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Exchange of L-Asparagine\n mean flux, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.000495, q= 0.0406")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Mean flux,\nmmol/(gDW * h)", limits=c(-80,0)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Exchange of L-Asparagine, sum
ggpaired(subset(sum_analyzed_fluxpaths_paired_wide, fluxpath_name %in% c("EX_asn_L(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Exchange of L-Asparagine\n sum of fluxes, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.000827, q= 0.0678")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Sum of fluxes,\nmmol/(gDW * h)", limits=c(-4000,0)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")





#############OLD STUFF BELOW####################


#########NEW Step 2: drop all nan values
analyzed_fluxpaths_noNANs <- analyzed_fluxpaths[!is.na(analyzed_fluxpaths$flux_value),]
#how many unique fluxes are we including here? 
fluxpaths_unique_noNANs <- unique(analyzed_fluxpaths_noNANs$fluxpath_name) #does not change # unique fluxpaths

#########Step 3: filter to only include paired taxa (i.e. taxa which appear in both tumor and normal samples)
analyzed_fluxpaths_wide <- select(analyzed_fluxpaths_noNANs, -c(Sample_ID)) #drop unnecessary columns
#first, convert to wide form
analyzed_fluxpaths_wide <- analyzed_fluxpaths_wide %>% 
  spread(Description, flux_value)

#then, drop all rows where both tumor and normal values are nan.
analyzed_fluxpaths_wide_noNANs <- analyzed_fluxpaths_wide[!(is.na(analyzed_fluxpaths_wide$tumor) & is.na(analyzed_fluxpaths_wide$normal)),]
fluxpaths_unique_zerofree <- unique(analyzed_fluxpaths_wide_noNANs$fluxpath_name) #count the number of metabolites
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
fluxpaths_zerofree_unique <- unique(analyzed_fluxpaths_wide_noNANs_zerofree$fluxpath_name)

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
#now, remove all rows where n (the name of the count column) is less than ...
analyzed_fluxpaths_filtered <- filter(analyzed_fluxpaths_filter_added, n >= 88)
#count number of unique metabolites left
fluxpaths_filtered_unique <- unique(analyzed_fluxpaths_filtered$fluxpath_subsystem)













#OOOOOOOLD STUFF
# let's look at SCFAs in particular: EX_but(e), EX_ac(e), 	EX_ppa(e), "EX_isobut(e)", "EX_4hphac(e)", "EX_acac(e)", "EX_pac(e)":
#butyrate, acetate, propionate, isobutyrate, 4-Hydroxyphenylacetate, Acetoacetate, Phenylacetate
scfa_names <- c("EX_but(e)","EX_ac(e)","EX_ppa(e)", "EX_isobut(e)", "EX_4hphac(e)", "EX_acac(e)", "EX_pac(e)")
analyzed_fluxpaths_scfas <- filter(analyzed_fluxpaths_noNANs, fluxpath_name %in% scfa_names) 

#How much does this change unique fluxpaths?
fluxpaths_unique_scfas <- unique(analyzed_fluxpaths_scfas$fluxpath_name)


#now, let's just... graph. y axis is flux, x-axis is patient_blind_id, tumor/normal samples graphed together.
ggplot(subset(analyzed_fluxpaths_scfas, fluxpath_name %in% c("EX_but(e)")), aes(x= Patient_Blind_ID, y= flux_value), fill = "Description") +
  geom_boxplot(aes(fill = Description), position = position_dodge(0.9)) +
  labs(title = expression("Butyrate exchange flux, Burns2015 data")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Flux,\nmmol/(gDW * h)", limits=c(-1000, 1000)) +
  #scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#this is fine, but I think we probably actually want to graph differences (more specifically, the per-taxon differences between tumor and normal samples)
#first, convert to wide form
analyzed_fluxpaths_scfas_wide <- select(analyzed_fluxpaths_scfas, -c(Sample_ID)) #drop unnecessary columns
analyzed_fluxpaths_scfas_wide <- analyzed_fluxpaths_scfas_wide %>% 
  spread(Description, flux_value)
#nans have been introduced- presumably because there are cases where a specific otu_id/flux combo is absent.
#that's probably fine. Add a difference column.
analyzed_fluxpaths_scfas_wide$difference <- (analyzed_fluxpaths_scfas_wide$tumor - analyzed_fluxpaths_scfas_wide$normal)
#now, remove all cols where difference=nan: 
analyzed_fluxpaths_scfas_wide_noNANs <- analyzed_fluxpaths_scfas_wide[!is.na(analyzed_fluxpaths_scfas_wide$difference),]
#go from 10434 rows to... 4348. Interesting. 

#let's try graphing again.
ggplot(subset(analyzed_fluxpaths_scfas_wide_noNANs, fluxpath_name %in% c("EX_pac(e)")), aes(x= Patient_Blind_ID, y= difference), fill = "Patient_Blind_ID") +
  geom_boxplot(aes(fill = Patient_Blind_ID), position = position_dodge(0.9)) +
  labs(title = expression("Butyrate exchange flux, Burns2015 data")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Difference in flux,\nmmol/(gDW * h)", limits=c(-1000, 1000)) +
  #scale_x_discrete(labels = c('Normal','Tumor')) +
  #scale_fill_brewer(palette = "Dark2") +
  #scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#it's colourful, but not actually that helpful. 

#all right, let's go with option 1. See if median vs mean has any effect. Then just look at specific metabolites maybe?
#EX_but(e), EX_ac(e), 	EX_ppa(e), "EX_isobut(e)", "EX_4hphac(e)", "EX_acac(e)", "EX_pac(e)":
analyzed_fluxpaths_noNANs
#add a sorting col
analyzed_fluxpaths_noNANs <- analyzed_fluxpaths_noNANs %>% 
  mutate(sorting_col = paste0(Sample_ID, "_", fluxpath_name)) 
#Take median, average, and summed fluxpath value across taxa, within a community
analyzed_fluxpaths_collapsed <- analyzed_fluxpaths_noNANs %>%        # Specify data frame
  group_by(sorting_col, Sample_ID, Description, Patient_Blind_ID, fluxpath_description, fluxpath_name) %>%    # Specify group indicator (columns to keep)
  summarise_at(vars(flux_value),                   # Specify column
               list(med_fluxpath_amount = median, av_fluxpath_amount=mean, sum_fluxpath_amount = sum)) %>% # Specify function
  ungroup() #do this so select (dropping columns) can happen later
#nice. End up with a dataframe with lots of graphing options for subsystems.
#let's do a graph.
ggplot(subset(analyzed_fluxpaths_collapsed, fluxpath_name %in% c("EX_pac(e)")), aes(x= Patient_Blind_ID, y= med_fluxpath_amount), fill = "Patient_Blind_ID") +
  geom_boxplot(aes(fill = Patient_Blind_ID), position = position_dodge(0.9)) +
  labs(title = expression("Butyrate exchange flux, Burns2015 data")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Difference in flux,\nmmol/(gDW * h)", limits=c(-1000, 1000)) +
  #scale_x_discrete(labels = c('Normal','Tumor')) +
  #scale_fill_brewer(palette = "Dark2") +
  #scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")
#terrible.

#let's plot median vs mean, tumor v normal
#ugh, need to convert to wide. double ugh, need to create separate dfs for each. 
analyzed_fluxpaths_med_wide <- select(analyzed_fluxpaths_collapsed, -c(sorting_col, Sample_ID, av_fluxpath_amount, sum_fluxpath_amount)) #drop unnecessary columns
analyzed_fluxpaths_med_wide <- analyzed_fluxpaths_med_wide %>% 
  spread(Description, med_fluxpath_amount)
analyzed_fluxpaths_med_wide$diff <- (analyzed_fluxpaths_med_wide$tumor - analyzed_fluxpaths_med_wide$normal)
analyzed_fluxpaths_med_wide <- analyzed_fluxpaths_med_wide %>%
  mutate(Difference_direction = case_when(
    diff > 0 ~ "tumor > normal",
    diff < 0 ~ "tumor < normal",
    diff == 0 ~ "Zero change"))
###
ggplot(subset(analyzed_fluxpaths_collapsed, fluxpath_name %in% c("EX_pac(e)")), aes(x= Description, y= med_fluxpath_amount), fill = "Description") +
  #geom_boxplot(aes(fill = Patient_Blind_ID), position = position_dodge(0.9)) +
  labs(title = expression("Butyrate exchange flux, Burns2015 data")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Difference in flux,\nmmol/(gDW * h)", limits=c(-1000, 1000)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")
####
analyzed_fluxpaths_med_wide <- select(analyzed_fluxpaths_collapsed, -c(sorting_col, Sample_ID, av_fluxpath_amount, sum_fluxpath_amount)) #drop unnecessary columns
analyzed_fluxpaths_med_wide <- analyzed_fluxpaths_med_wide %>% 
  spread(Description, med_fluxpath_amount)
analyzed_fluxpaths_med_wide$diff <- (analyzed_fluxpaths_med_wide$tumor - analyzed_fluxpaths_med_wide$normal)
analyzed_fluxpaths_med_wide <- analyzed_fluxpaths_med_wide %>%
  mutate(Difference_direction = case_when(
    diff > 0 ~ "tumor > normal",
    diff < 0 ~ "tumor < normal",
    diff == 0 ~ "Zero change"))
###
ggpaired(subset(analyzed_fluxpaths_med_wide, fluxpath_name %in% c("EX_pac(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Fatty acid oxidation subsystem\n median flux, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.00599, q= 0.159")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Median flux,\nmmol/(gDW * h)", limits=c(-1000,1000)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")
#######all right, seems like median ain't it. Too many median values are 0.
analyzed_fluxpaths_av_wide <- select(analyzed_fluxpaths_collapsed, -c(sorting_col, Sample_ID, med_fluxpath_amount, sum_fluxpath_amount)) #drop unnecessary columns
analyzed_fluxpaths_av_wide <- analyzed_fluxpaths_av_wide %>% 
  spread(Description, av_fluxpath_amount)
analyzed_fluxpaths_av_wide$diff <- (analyzed_fluxpaths_av_wide$tumor - analyzed_fluxpaths_av_wide$normal)
analyzed_fluxpaths_av_wide <- analyzed_fluxpaths_av_wide %>%
  mutate(Difference_direction = case_when(
    diff > 0 ~ "tumor > normal",
    diff < 0 ~ "tumor < normal",
    diff == 0 ~ "Zero change"))
###
ggpaired(subset(analyzed_fluxpaths_av_wide, fluxpath_name %in% c("EX_ppa(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Fatty acid oxidation subsystem\n median flux, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.00599, q= 0.159")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Median flux,\nmmol/(gDW * h)", limits=c(0,1000)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")
#######Seems like average isn't it either. differences are too small.
analyzed_fluxpaths_sum_wide <- select(analyzed_fluxpaths_collapsed, -c(sorting_col, Sample_ID, med_fluxpath_amount, av_fluxpath_amount)) #drop unnecessary columns
analyzed_fluxpaths_sum_wide <- analyzed_fluxpaths_sum_wide %>% 
  spread(Description, sum_fluxpath_amount)
analyzed_fluxpaths_sum_wide$diff <- (analyzed_fluxpaths_sum_wide$tumor - analyzed_fluxpaths_sum_wide$normal)
analyzed_fluxpaths_sum_wide <- analyzed_fluxpaths_sum_wide %>%
  mutate(Difference_direction = case_when(
    diff > 0 ~ "tumor > normal",
    diff < 0 ~ "tumor < normal",
    diff == 0 ~ "Zero change"))
###
ggpaired(subset(analyzed_fluxpaths_sum_wide, fluxpath_name %in% c("EX_ppa(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Fatty acid oxidation subsystem\n median flux, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.00599, q= 0.159")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Median flux,\nmmol/(gDW * h)", limits=c(-0,30000)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")
#EX_but(e), EX_ac(e), 	EX_ppa(e), "EX_isobut(e)", "EX_4hphac(e)", "EX_acac(e)", "EX_pac(e)":

#similarly uninformative

#what if we plotted PCA instead?
#go back to analyzed_fluxpaths_noNANs and convert to wide...
analyzed_fluxpaths_noNANs_wide <- select(analyzed_fluxpaths_noNANs, -c(sorting_col, fluxpath_description)) #drop unnecessary columns
analyzed_fluxpaths_noNANs_wide <- analyzed_fluxpaths_noNANs_wide %>% 
  spread(fluxpath_name, flux_value)

#allrighty, that's that. load our libraries for pca
library(tidyverse)
library(broom)

#NOTE: we need to set a few variables as categorical here.
analyzed_fluxpaths_noNANs_wide$OTU_ID = as.factor(analyzed_fluxpaths_noNANs_wide$OTU_ID)
analyzed_fluxpaths_noNANs_wide$Sample_ID = as.factor(analyzed_fluxpaths_noNANs_wide$Sample_ID)
analyzed_fluxpaths_noNANs_wide$Patient_Blind_ID = as.factor(analyzed_fluxpaths_noNANs_wide$Patient_Blind_ID)

#"We use select() to select numerical variables in penguins’s data, apply scale() and then do PCA with prcomp() function."
pca_fit <- analyzed_fluxpaths_noNANs_wide %>%
  select(where(is.numeric)) %>%
  scale() %>%
  prcomp()

#will NOT WORK because we have missing data. Cool.See http://juliejosse.com/wp-content/uploads/2018/05/DataAnalysisMissingR.html
analyzed_fluxpaths_noNANs_wide_quantonly <- select (analyzed_fluxpaths_noNANs_wide, -c(OTU_ID, Sample_ID, Patient_Blind_ID, Description))

nb <- estim_ncpPCA(analyzed_fluxpaths_noNANs_wide_quantonly,method.cv = "Kfold", verbose = FALSE)
#NOTE: this takes forever. 

nb$ncp
res.comp <- imputePCA(don, ncp = nb$ncp) # iterativePCA algorithm
res.comp$completeObs[1:3,] # the imputed data set

#don't know if this will work, but anyway.
pca_fit <- res.comp %>%
  #select(where(is.numeric)) %>%
  scale() %>%
  prcomp()


#Here is a quick summary of the PCA. We can see that the first two principal components explain 88% of the variation in the data.
summary(pca_fit)



