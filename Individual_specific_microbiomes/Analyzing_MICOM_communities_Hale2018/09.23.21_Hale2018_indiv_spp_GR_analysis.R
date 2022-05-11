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

#further based on Date_TBA_microbial_GR_analysis_Burns_data_filtered_top_20_taxa.R, with updates for Hale data
#and 09.03.21_Hale2018_indiv_spp_GR_analysis.R

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Hale2018/") 

##########Step 1: import data and adjust factors
#It is important to note here that I have data (and have done analyses) for both family and genus level. However, I went with genus level
#moving forward because... most of the family level associations are basically due to a single genus anyway, I think is what I decided. 
#NOTE: this is a temporary file; edit as needed
GR_by_genus_point_one <- read_csv("09.01.21_Hale2018_indiv_spp_GRs_collapsed_to_genus_diffs.csv")
GR_by_genus_point_one$host_subject_id <- as.factor(GR_by_genus_point_one$host_subject_id)
#colnames(GR_by_genus_point_one) #print column names as needed

#########Step 2: reorder and add columns for filtering and stats
#add a new column, famgen, for sorting purposes- some genera have the same name but different families. 
GR_by_genus_point_one <- GR_by_genus_point_one  %>% unite("famgen_col", Family:Genus, remove = FALSE) 
colnames(GR_by_genus_point_one) #check
#now rearrange so famgen and host_subject_id are next to each other. 
GR_by_genus_point_one <- GR_by_genus_point_one[, c(1,2,3,4,5,6,19,7,8,9,10,11,12,13,14,15,16,17,18,20,21,22,23,24)] #this is edited from previous version, as my output is a bit different

#how many unique taxa are we including here? As of 09.22.21, 130 taxa
GR_by_genus_unique <- unique(GR_by_genus_point_one$famgen_col)

#########Step 3 NEW: filter by taxa which appear in at least 20% of the total number of samples 
#again, borrow from Sambhawa's code for filtering genes; see 05.03.21_DEseq2_analysis_Burns_data_filtered.R
#Filter to only include taxa which appear in more than a quarter of the samples, which translates to... number of times a value appears in famgen has to be 
#greater than or equal to 25 
#NEW in this version: remove all taxa samples (in this dataframe, rows) where mean_genus_GR_normal and mean_genus_GR_tumor are both zero
#fortunately, we have a built in way to do this: any rows where mean_GR_difference is 0 should be removed. 
GR_by_genus_point_one_zerofree <- GR_by_genus_point_one[(GR_by_genus_point_one$mean_GR_difference !=0),]
#all right, by this point, we know that each row represents paired samples; by definition, each taxon must appear in at least two samples
#now, we need some way to count the number of samples where each famgen col appears. 
#this is different than counting the number of occurances of each famgen value; 
#uhh, scratch that, no it isn't, because we've collapsed the data already such that each row represents the mean GR of that taxon in that sample
#meaning if there were multiple instances of that taxon within a sample, it's already been averaged out. 
#so, our method is sound. count the number of each occurance of famgen. Multiply that by 2. That is the number of SAMPLES (not patients)
#in which that taxon appears. 

#how many unique taxa are we including here? As of 09.22.21, 130 taxa
GR_by_genus_zerofree_unique <- unique(GR_by_genus_point_one_zerofree$famgen_col)

#then, create a filtering dataframe that contains counts of each famgen incidence in the dataframe
GR_by_genus_point_one_filtering_df <- GR_by_genus_point_one_zerofree %>%  #create a new dataframe by filtering the old one
  mutate(count_of_nonzeros = rowSums(.[10:11]!=0)) %>% #create a new column, count_of_nonzeroes, that counts # nonzero values in columns 10 and 11 for each row
  group_by(famgen_col) %>% #group by the column of interest
  summarize(n = sum(count_of_nonzeros)) %>% #sums the count_of_nonzeros values within each famgen and saves the value as n
  mutate(n) #add n as a new column in the dataframe

#create a new dataframe that merges the old and new dataframes
GR_by_genus_point_one_filter_added <- dplyr::inner_join(GR_by_genus_point_one_zerofree, GR_by_genus_point_one_filtering_df, by= "famgen_col")

#now, remove all rows where n (the name of the count column) is less than 18 (20% of the total number of samples)
GR_by_genus_point_one_filtered <- filter(GR_by_genus_point_one_filter_added, n >= 18)
GR_by_genus_filtered_unique <- unique(GR_by_genus_point_one_filtered$famgen_col) #count unique taxa
GR_by_genus_filtered_unique

rm(GR_by_genus_point_one, GR_by_genus_point_one_filter_added, GR_by_genus_point_one_zerofree, GR_by_genus_point_one_filtering_df)

#########Step 4: filter to only include paired taxa (i.e. taxa which appear in both tumor and normal samples)
#create a new dataframe where only taxa with growth rates in BOTH tumor and normal samples (i.e. host_subject_id appears twice) are included
#genus_point_one_paired <- GR_by_genus_point_one_filtered   %>% unite("sorting_col", famgen_col:host_subject_id, remove = FALSE)
#genus_point_one_paired_only <- subset(genus_point_one_paired,duplicated(sorting_col) | duplicated(sorting_col, fromLast=TRUE))
#after all this filtering, we get 1568 taxa-host_subject_id combinations. Not sure how this compares to previous results
#I'm not sure why... we're doing this. By definition the data are already paired. 
#genus_point_one_paired_only <- copy(GR_by_genus_point_one_filtered)
#oh, I think there was a mix-up between wide and longform data here.

##########Step 5: split into multiple dataframes by genus and do stats on each genus
#first, edit the dataframe
#drop a couple of columns
genus_point_one_paired_only <- select(GR_by_genus_point_one_filtered, -c(...1, mean_GR_difference, n, Description))
#rename a couple of columns to make life easier later. Yeah, yeah, I know, this is the reverse of what I did before.
genus_point_one_paired_only <- genus_point_one_paired_only %>% rename(normal= mean_genus_GR_normal, tumor= mean_genus_GR_tumor) #rename these columns
#convert to long form
genus_point_one_paired_only <- reshape2::melt(data= genus_point_one_paired_only,
                                    id.vars= c("Kingdom", "Phylum", "Class", "Order", "famgen_col", "host_subject_id", "Family", "Genus",
                                    "SampleID", "cancer_stage", "env_material", "ethnicity", "Host_Age", "host_body_mass_index", "host_sex",
                                    "host_tissue_sampled", "mismatch_repair_status", "collection_site", "biopsy_number"),
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
rm(elem, genus_0.1_stats, genus_point_one_paired, genus_point_one_paired_only, split_0.1_GR_genus, 
   genus_0.1_pvals, genus_0.1_qvals, genus_ids_point_one, new_value)

###########Step 8: pull out significant p-values only
genus_0.1_statistics <- genus_0.1_statistics[order(genus_0.1_statistics$genus_0.1_pvals),] #re-order by p-value
#top_ten_genera <- genus_0.1_statistics %>% top_n(-20, genus_0.1_pvals) #note the -20 here to select the smallest p-values (20 taxa)
#EYEEEEE, there are only nine taxa here, but at least the q- values actually mean something!

##########Step 9: save results as a csv file
write.csv(as.data.frame(genus_0.1_statistics), file="09.23.21_filtered_microbial_GR_all_results_Hale2018_data.csv")
#new: filter by sig p-values and only keep significantly different p-values
genus_0.1_statistics_sig <- genus_0.1_statistics[genus_0.1_statistics$genus_0.1_pvals <= 0.05, ]
write.csv(as.data.frame(genus_0.1_statistics_sig), file="09.23.21_filtered_microbial_GR_SIG_results_Hale2018_data.csv")


#DID NOT DO GRAPHING HERE
##########Step 10: GRAPH??!!!
GR_by_genus_wide <- GR_by_genus_point_one_filtered %>%
  mutate(Difference_direction = case_when(
    mean_GR_difference > 0 ~ "tumor > normal",
    mean_GR_difference < 0 ~ "tumor < normal",
    mean_GR_difference == 0 ~ "Zero change",
  ))

write.csv(as.data.frame(GR_by_genus_wide), file="09.23.21_filtered_microbial_GR_data_all_Hale2018_data.csv")


#GRAPH!!!

#count number of samples with Fuso:
freqs <- table(GR_by_genus_wide$Genus)
freqs

#Fusobacterium
ggpaired(subset(GR_by_genus_wide, Genus %in% c("Fusobacterium")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Fusobacterium growth rates, Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.000363, q= 0.0105")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.6)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Bacteroides
ggpaired(subset(GR_by_genus_wide, Genus %in% c("Bacteroides")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Bacteroides growth rates, Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.00819, q= 0.119")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.6)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Negativibacillus
ggpaired(subset(GR_by_genus_wide, Genus %in% c("Negativibacillus")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Negativibacillus growth rates, Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.0144, q= 0.14019")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.08)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Blautia
ggpaired(subset(GR_by_genus_wide, Genus %in% c("Blautia")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Blautia growth rates, Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.0313, q= 0.222")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.3)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Subdoligranulum
ggpaired(subset(GR_by_genus_wide, Genus %in% c("Subdoligranulum")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Subdoligranulum growth rates, Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.0383, q= 0.222")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.25)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Parabacteroides
ggpaired(subset(GR_by_genus_wide, Genus %in% c("Parabacteroides")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Parabacteroides growth rates, Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.0489, q= 0.263")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.6)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Faecalibacterium
ggpaired(subset(GR_by_genus_wide, Genus %in% c("Faecalibacterium")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Faecalibacterium growth rates, Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.0858, q= 0.129")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.6)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Collinsella
ggpaired(subset(GR_by_genus_wide, Genus %in% c("Collinsella")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Collinsella growth rates, Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.107, q= 0.138")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.1)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Escherichia/Shigella
ggpaired(subset(GR_by_genus_wide, Genus %in% c("Escherichia/Shigella")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Escherichia/Shigella growth rates, Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.770, q= 0.866")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,1.6)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Agathobacter
ggpaired(subset(GR_by_genus_wide, Genus %in% c("Agathobacter")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Agathobacter growth rates, Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.926, q= 0.926")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.6)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")


