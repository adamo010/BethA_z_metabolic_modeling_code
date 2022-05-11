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

#the goal here is to calculate p-values for the significant differences between input/output metabolites to use for correlating with gene expression
#differences. Also, we want to abundance filter input/output metabolites that do not appear in at least half the number of samples 
#(before running stats on identifying metabolites with significantly different input/outputs between tumor and normal samples).

#loosely based on 05.14.21_microbial_GR_analysis_Burns_data_filtered_top_20_taxa.R

#further based on Date_TBA_microbial_inputoutput_analysis_Burns_data_filtered_top_20_metabolites.R, and 09.07.21_Hale2018_inputoutput_analysis.R
#same as 09.10.21_Burns2015_inputoutput_analysis.R, but with added step of collapsing trace minerals components.

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/") 

##########Step 1: import data and adjust factors
#It is important to note here that I have data (and have done analyses) for both family and genus level. However, I went with genus level
#moving forward because... most of the family level associations are basically due to a single genus anyway, I think is what I decided. 
analyzed_metabolites <- read_csv("09.09.21_Burns2015_inputoutput_with_metadata_longform.csv")
analyzed_metabolites$Patient_Blind_ID <- as.factor(analyzed_metabolites$Patient_Blind_ID)
#colnames(analyzed_metabolites) #print column names as needed

#########Step 2 NEW as of 09.20.21: collapsing by trace minerals
#there are several metabolites (Calcium, Chloride ion, Cobalt, Cu2+, Heme, Magnesium, Mn2+, Potassium, Siroheme, and Zinc (II) ion)
#which have identical fluxes in all communities. The goal is to collapse these into a single group, called "Trace minerals"
#first, create a vector of these metabolites
trace_minerals <- c("Calcium", "Chloride ion", "Cobalt", "Cu2+", "Heme", "Magnesium", "Mn2+", "Potassium", "Siroheme", "Zinc (II) ion")
#Figure out some way to count the number of trace minerals with nonzero values within each sample. 
#As long as it's always 10 or 0, can collapse. 
#what the hell, I got rid of my sample_id column. make a new one (of sorts)
analyzed_metabolites <- analyzed_metabolites %>% unite("sample_ID_ish", Description:Patient_Blind_ID, remove = FALSE)
#first, create a dataframe that ONLY contains trace minerals
analyzed_metabolites_trace_minerals_only <- analyzed_metabolites[analyzed_metabolites$metabolite_name %in% trace_minerals, ]

#then, create a filtering_DF that counts the number of unique values within each sample_ID_ish
trace_minerals_filtering_df <- analyzed_metabolites_trace_minerals_only %>%  #create a new dataframe by filtering the old one
  group_by(sample_ID_ish) %>% #group by the column of interest
  summarize(n = n_distinct(metabolite_amount)) #count the number of occurrences of each unique value in metabolite_amount and store as variable 'n'

#damn it. It's a range from 3-10.
#welll, is it? I took a look at 2343, which allegedly has three distinct amounts, but it sure looks like one amount.
#I bet this is a rounding issue. Use dplyr to round to 6 significant digits. That seems good enough
analyzed_metabolites_trace_minerals_only <- analyzed_metabolites_trace_minerals_only %>% 
  mutate(metabolite_amount=signif(metabolite_amount, 6)) 

#rerun the filtering_DF
trace_minerals_filtering_df <- analyzed_metabolites_trace_minerals_only %>%  #create a new dataframe by filtering the old one
  group_by(sample_ID_ish) %>% #group by the column of interest
  summarize(n = n_distinct(metabolite_amount)) #count the number of occurrences of each unique value in metabolite_amount and store as variable 'n'

#okay, now we have all 1s, with a few 2s. This is actually good. It means within each sample_ID_ish, there is only 1-2 unique
#metabolite_amounts. Those which have 2 unique amounts have ###, and 0. 
#it seems like there are 21 samples which have zero values. These are zeroes for heme and seroheme.So... let's just ignore it.

#time to collapse! for each sample_ID_ish, only keep ONE of the rows where metabolite_names values are in the trace_minerals vector
#and rename metabolite_name to trace_mineral.
analyzed_metabolites_collapsed <-data.frame(analyzed_metabolites) #copy the original dataframe
#first, rename "Potassium" to "trace elements"
analyzed_metabolites_collapsed$metabolite_name[analyzed_metabolites_collapsed$metabolite_name== "Potassium"] <- "Trace_elements"
#then, remove all rows where trace_minerals vector values are in the metabolite_name column
analyzed_metabolites_collapsed<- analyzed_metabolites_collapsed[!(analyzed_metabolites_collapsed$metabolite_name %in% trace_minerals),]

#NICE! Let's move forward. 

#########Step 2: filter by metabolites which appear in at least half the total number of samples 
#we'll set a cutoff for at least 10-6.
#Filter to only include metabolites which appear in more than half of the samples, which translates to... 
#number of times a value appears greater than 10-6 in metabolite_name has to be greater than or equal to 44
#first, create a filtering dataframe that contains counts of each metabolite_name incidence in the dataframe
analyzed_metabolites_filtering_df <- analyzed_metabolites_collapsed %>%  #create a new dataframe by filtering the old one
  group_by(metabolite_name) %>% #group by the column of interest
  summarize(n = n_distinct(metabolite_amount)) #count the number of occurrences of each unique value in metabolite_amount and store as variable 'n'
#subtracting 1 b/c nan is counted as a unique value
#NOTE: use this dataframe to see how many taxa would be included at each cutoff

#create a new dataframe that merges the old and new dataframes
analyzed_metabolites_filter_added <- dplyr::inner_join(analyzed_metabolites, analyzed_metabolites_filtering_df, by= "metabolite_name")
#now, remove all rows where n (the name of the count column) is less than 13
analyzed_metabolites_filtered <- filter(analyzed_metabolites_filter_added, n >= 13)
rm(analyzed_metabolites_filtering_df, analyzed_metabolites_filter_added)

#########Step 3: filter to only include paired taxa (i.e. taxa which appear in both tumor and normal samples)
#first, drop all rows that contain nan in metabolite_amount column
analyzed_metabolites_filtered_noNANs <- analyzed_metabolites_filtered[!is.na(analyzed_metabolites_filtered$metabolite_amount),]
#then, do some reordering- maybe this is just a Hale data thing? check if Patient_Blind_id is next to metabolite_name
analyzed_metabolites_filtered_noNANs <- analyzed_metabolites_filtered_noNANs[, c(1,2,6,3,4,5)] #this is edited from previous version, as my output is a bit different
#then, create a new dataframe where only taxa with growth rates in BOTH tumor and normal samples (i.e. Patient_Blind_id appears twice) are included
analyzed_metabolites_paired <- analyzed_metabolites_filtered_noNANs %>% unite("sorting_col", metabolite_name:Patient_Blind_ID, remove = FALSE)
analyzed_metabolites_paired_only <- subset(analyzed_metabolites_paired,duplicated(sorting_col) | duplicated(sorting_col, fromLast=TRUE))
#after all this filtering, we get 1006 metabolite-Patient_blind_id combinations. 

##########Step 4: split into multiple dataframes by metabolite and do stats on each metabolite
split_analyzed_metabolites_by_metab <- split(analyzed_metabolites_paired_only, with(analyzed_metabolites_paired_only, interaction(metabolite_name)), drop = TRUE)
by_metab_stats <- lapply(split_analyzed_metabolites_by_metab, function(df){wilcox.test(metabolite_amount~Description, data=df, exact= FALSE, paired= TRUE)})

##########Step 5: statistics
#make lists of p values
by_metab_pvals <- c() #initialize list
for (elem in by_metab_stats){
  new_value = elem$p.value
  by_metab_pvals <- c(by_metab_pvals, new_value)}
rm(new_value)
#make a list of q-values
by_metab_qvals <- p.adjust(by_metab_pvals, method = "fdr")
#make a list of genus_ids
metab_ids <- c()
for(elem in split_analyzed_metabolites_by_metab){
  new_value = elem$metabolite_name[1]
  metab_ids <- c(metab_ids, new_value)}
#merge all lists together. 
metabolite_statistics <- data.frame(metab_ids, by_metab_pvals, by_metab_qvals)

###########Step 6: clean up and remove extra variables
rm(elem, analyzed_metabolites_paired, split_analyzed_metabolites_by_metab, 
   by_metab_pvals, by_metab_qvals, metab_ids, new_value)

###########Step 7: pull out significant p-values only
metabolite_statistics <- metabolite_statistics[order(metabolite_statistics$by_metab_pvals),] #re-order by p-value
#keep all 22. Whatever. 

##########Step 8: save results as a csv file
write.csv(as.data.frame(metabolite_statistics), file="09.21.21_filtered_microbial_input_output_sig_results_Burns2015_data_trace_elements_combined.csv")

##########Step 9: prepare data for graphing START HERE
#LAAAAAAAAAAAAME now I have to calculate differences.
#use analyzed_metabolites_paired_only. First, drop boring cols
analyzed_metabolites_paired_wide = subset(analyzed_metabolites_paired_only, select = -c(...1, sample_ID_ish, sorting_col) )
#Split up to wide.
analyzed_metabolites_paired_wide <- spread(analyzed_metabolites_paired_wide, Description, metabolite_amount)
#add a new column called diff. tumor-normal
analyzed_metabolites_paired_wide$diff <- (analyzed_metabolites_paired_wide$tumor - analyzed_metabolites_paired_wide$normal)

#NOW we can add a new column with useful binning of differences for graphing purposes
analyzed_metabolites_paired_wide <- analyzed_metabolites_paired_wide %>%
  mutate(Difference_direction = case_when(
    diff > 0 ~ "tumor > normal",
    diff < 0 ~ "tumor < normal",
    diff == 0 ~ "Zero change",
  ))

write.csv(as.data.frame(analyzed_metabolites_paired_wide), file="09.21.21_filtered_microbial_input_output_data_Burns2015_data_trace_elements_combined.csv")

#step 10: NOW GRAPH??!!!
#Adenosine monophosphate
ggpaired(subset(analyzed_metabolites_paired_wide, metabolite_name %in% c("Adenosine monophosphate")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Adenosine monophosphate import, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.0176, q= 0.105")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(0,0.01)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Fe3+
ggpaired(subset(analyzed_metabolites_paired_wide, metabolite_name %in% c("Fe3+")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Fe3+ import, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.132, q= 0.397")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(0,0.00012)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#D-Galactose
ggpaired(subset(analyzed_metabolites_paired_wide, metabolite_name %in% c("D-Galactose")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("D-Galactose export, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.281, q= 0.561")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(-0.006,0)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#mucin-type O-glycan No 9
ggpaired(subset(analyzed_metabolites_paired_wide, metabolite_name %in% c("mucin-type O-glycan No 9")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("mucin-type O-glycan No 9 import, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.554, q= 0.668")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(0, 0.006)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#released mucin-type O-glycan No 9
ggpaired(subset(analyzed_metabolites_paired_wide, metabolite_name %in% c("released mucin-type O-glycan No 9")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("released mucin-type O-glycan No 9 import, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.610, q= 0.668")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(0, 0.006)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Fe2+
ggpaired(subset(analyzed_metabolites_paired_wide, metabolite_name %in% c("Fe2+")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Fe2+ import, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.668, q= 0.668")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(0, 0.00012)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")
