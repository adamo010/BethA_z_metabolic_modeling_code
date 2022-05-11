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

#the goal here is to calculate p-values for the significant differences between taxa growth rates to use for correlating with gene expression
#differences. Also, we want to abundance filter taxa that do not appear in at least half the number of samples (before running stats on 
#identifying taxa with significantly different growth rates between tumor and normal samples).

#loosely based on 04.06.21_correlating_RNAseq_and_GRs_Burns_data.R, but without the correlating part; that goes in the next script
#an updated version of 05.14.21_microbial_GR_analysis_Burns_data_filtered_top_20_taxa.R, with new file paths
#and new data from rerunning MICOM on MSI with an updated version.

#further based on Date_TBA_microbial_GR_analysis_Burns_data_filtered_top_20_taxa.R and 09.03.21_Hale2018_indiv_spp_GR_analysis.R

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/") 

##########Step 1: import data and adjust factors
#It is important to note here that I have data (and have done analyses) for both family and genus level. However, I went with genus level
#moving forward because... most of the family level associations are basically due to a single genus anyway, I think is what I decided. 
#NOTE: this is a temporary file; edit as needed
GR_by_genus_point_one <- read_csv("09.09.21_Burns2015_indiv_spp_GRs_collapsed_to_genus_diffs.csv")
GR_by_genus_point_one$Patient_Blind_ID <- as.factor(GR_by_genus_point_one$Patient_Blind_ID)
#colnames(GR_by_genus_point_one) #print column names as needed

#########Step 2: reorder and add columns for filtering and stats
#add a new column, famgen, for sorting purposes- some genera have the same name but different families. 
GR_by_genus_point_one <- GR_by_genus_point_one  %>% unite("famgen_col", family:genus, remove = FALSE) 
colnames(GR_by_genus_point_one) #check

#now rearrange so famgen and Patient_Blind_ID are next to each other. 
GR_by_genus_point_one <- GR_by_genus_point_one[, c(1,2,3,4,5,6,14,7,8,9,10,11,12,13,15)] #this is edited from previous version, as my output is a bit different

#########Step 3 NEW: filter by taxa which appear in at least 20% of the total number of samples 
#again, borrow from Sambhawa's code for filtering genes; see 05.03.21_DEseq2_analysis_Burns_data_filtered.R
#Filter to only include taxa which appear in more than a quarter of the samples, which translates to... number of times a value appears in famgen has to be 
#greater than or equal to 25 
#NEW in this version: remove all taxa samples (in this dataframe, rows) where mean_genus_GR_normal and mean_genus_GR_tumor are both zero
#fortunately, we have a built in way to do this: any rows where mean_GR_difference is 0 should be removed. 
GR_by_genus_point_one_zerofree <- GR_by_genus_point_one[(GR_by_genus_point_one$mean_GR_difference !=0),]

#then, create a filtering dataframe that contains counts of each famgen incidence in the dataframe
GR_by_genus_point_one_filtering_df <- GR_by_genus_point_one_zerofree %>%  #create a new dataframe by filtering the old one
  group_by(famgen_col) %>% #group by the column of interest
  summarize(n = n()) %>% #count the number of occurrences of each unique value in famgen_col and store as variable 'n'
  mutate(n) #add n as a new column in the dataframe
#NOTE: use this dataframe to see how many taxa would be included at each cutoff

#create a new dataframe that merges the old and new dataframes
GR_by_genus_point_one_filter_added <- dplyr::inner_join(GR_by_genus_point_one, GR_by_genus_point_one_filtering_df, by= "famgen_col")
#now, remove all rows where n (the name of the count column) is less than 11
GR_by_genus_point_one_filtered <- filter(GR_by_genus_point_one_filter_added, n >= 13)
rm(GR_by_genus_point_one_filtering_df, GR_by_genus_point_one_filter_added)

##########Step 5: split into multiple dataframes by genus and do stats on each genus
#first, edit the dataframe
#drop a couple of columns
genus_point_one_paired_only <- select(GR_by_genus_point_one_filtered, -c(...1, mean_GR_difference, n, Description))
#rename a couple of columns to make life easier later. Yeah, yeah, I know, this is the reverse of what I did before.
genus_point_one_paired_only <- genus_point_one_paired_only %>% rename(normal= mean_genus_GR_normal, tumor= mean_genus_GR_tumor) #rename these columns
#convert to long form
colnames(genus_point_one_paired_only) #check column names for melting
genus_point_one_paired_only <- reshape2::melt(data= genus_point_one_paired_only,
                                              id.vars= c("kingdom", "phylum", "class", "order", "famgen_col", "Patient_Blind_ID", "family", "genus",
                                                         "SampleID", "Site"),
                                              variable.name = "Description",
                                              value.name = "mean_genus_GR")
#split up by different genera and see where we get. 
split_0.1_GR_genus <- split(genus_point_one_paired_only, with(genus_point_one_paired_only, interaction(famgen_col)), drop = TRUE)
genus_0.1_stats <- lapply(split_0.1_GR_genus, function(df){wilcox.test(mean_genus_GR~Description, data=df, exact= FALSE, paired= TRUE)})

##########Step 6: statistics
#make lists of p values
genus_0.1_pvals <- c() #initialize list
for (elem in genus_0.1_stats){
  new_value = elem$p.value
  genus_0.1_pvals <- c(genus_0.1_pvals, new_value)}
rm(new_value)
#make a list of q-values
genus_0.1_qvals <- p.adjust(genus_0.1_pvals, method = "fdr")
#make a list of genus_ids
genus_ids_point_one <- c()
for(elem in split_0.1_GR_genus){
  new_value = elem$famgen_col[1]
  genus_ids_point_one <- c(genus_ids_point_one, new_value)}
#merge all lists together. 
genus_0.1_statistics <- data.frame(genus_ids_point_one, genus_0.1_pvals, genus_0.1_qvals)

###########Step 7: clean up and remove extra variables
rm(elem, genus_0.1_stats, genus_point_one_paired_only, split_0.1_GR_genus, 
   genus_0.1_pvals, genus_0.1_qvals, genus_ids_point_one, new_value)

###########Step 8: pull out significant p-values only
genus_0.1_statistics <- genus_0.1_statistics[order(genus_0.1_statistics$genus_0.1_pvals),] #re-order by p-value
#top_ten_genera <- genus_0.1_statistics %>% top_n(-20, genus_0.1_pvals) #note the -20 here to select the smallest p-values (20 taxa)
#hey, 20 samples exactly. Cool beans.

##########Step 9: save results as a csv file
write.csv(as.data.frame(genus_0.1_statistics), file="09.10.21_filtered_microbial_GR_sig_results_Burns2015_data.csv")

##########Step 10: GRAPH??!!!
GR_by_genus_wide <- GR_by_genus_point_one_filtered %>%
  mutate(Difference_direction = case_when(
    mean_GR_difference > 0 ~ "tumor > normal",
    mean_GR_difference < 0 ~ "tumor < normal",
    mean_GR_difference == 0 ~ "Zero change",
  ))

write.csv(as.data.frame(GR_by_genus_wide), file="09.10.21_filtered_microbial_GR_data_Burns2015_data.csv")

#GRAPH!!!
#Fusobacterium
ggpaired(subset(GR_by_genus_wide, genus %in% c("Fusobacterium")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Fusobacterium growth rates, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.00162, q= 0.0324")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.25)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Ruminococcus
ggpaired(subset(GR_by_genus_wide, genus %in% c("Ruminococcus")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Ruminococcus growth rates, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.0131, q= 0.108")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.3)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Lactobacillus
ggpaired(subset(GR_by_genus_wide, genus %in% c("Lactobacillus")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Lactobacillus growth rates, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.0163, q= 0.108")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.15)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Ralstonia
ggpaired(subset(GR_by_genus_wide, genus %in% c("Ralstonia")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Ralstonia growth rates, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.0832, q= 0.317")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.04)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Clostridium
ggpaired(subset(GR_by_genus_wide, genus %in% c("Clostridium")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Clostridium growth rates, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.0894, q= 0.317")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.4)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Faecalibacterium
ggpaired(subset(GR_by_genus_wide, genus %in% c("Faecalibacterium")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Faecalibacterium growth rates, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.0950, q= 0.317")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.6)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")




