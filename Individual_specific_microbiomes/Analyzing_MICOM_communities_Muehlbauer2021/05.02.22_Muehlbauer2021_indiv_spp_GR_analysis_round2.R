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

#the goal here is to... filter taxa, and see if anything interesting comes from the growth rates. NOT doing tumor v normal, as these don't apply here

#based on 05.02.22_Muehlbauer2021_indiv_spp_GR_analysis_round2.R

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Muehlbauer2021/") 

##########Step 1: import data and adjust factors
#It is important to note here that I have data (and have done analyses) for both family and genus level. However, I went with genus level
#moving forward because... most of the family level associations are basically due to a single genus anyway, I think is what I decided. 
#NOTE: this is a temporary file; edit as needed
GR_by_genus_point_one <- read_csv("05.02.22_Muehlbauer2021_indiv_spp_GRs_collapsed_to_genus.csv")
GR_by_genus_point_one <- GR_by_genus_point_one %>% select(-c(...1))
GR_by_genus_point_one$Isolate <- as.factor(GR_by_genus_point_one$Isolate)
#colnames(GR_by_genus_point_one) #print column names as needed

#ALL RIGHT, I'm making a game time decision here. We're ONLY keeping Colonocyte info. This will change downstream analysis, but I think we can
#handle it. 
GR_by_genus_point_one_colonly <- GR_by_genus_point_one[(GR_by_genus_point_one$experimental_treatment=="Colonocytes"),]
GR_by_genus_point_one_colonly <- GR_by_genus_point_one_colonly %>% select(-c(experimental_treatment, Library_name))

#clean up
rm(GR_by_genus_point_one)

#########Step 2: reorder and add columns for filtering and stats
#add a new column, famgen, for sorting purposes- some genera have the same name but different families. 
GR_by_genus_point_one_colonly <- GR_by_genus_point_one_colonly  %>% unite("famgen_col", family:genus, remove = FALSE) 
colnames(GR_by_genus_point_one_colonly) #check
#now rearrange so famgen and Isolate are next to each other. 
GR_by_genus_point_one_colonly <- GR_by_genus_point_one_colonly[, c(1,2,3,4,5,6,10,7,8,9)] #this is edited from previous version, as my output is a bit different

#how many unique taxa are we including here? As of 05.02.22, have 82 taxa
GR_by_genus_unique <- unique(GR_by_genus_point_one_colonly$famgen_col)

#########Step 3 NEW: filter by taxa which appear in at least 2 patients 
#remove all taxa where mean_genus_GR is 0
GR_by_genus_point_one_zerofree <- GR_by_genus_point_one_colonly[(GR_by_genus_point_one_colonly$mean_genus_GR !=0),]

#how many unique taxa are we including here? As of 05.02.22, 74 taxa
GR_by_genus_zerofree_unique <- unique(GR_by_genus_point_one_zerofree$famgen_col)

colnames(GR_by_genus_point_one_zerofree) #check mean_genus_GR_normal and mean_genus_GR_tumor column numbers; will be used below
#here, 9 and 10

#then, create a filtering dataframe that contains counts of each famgen incidence in the dataframe
GR_by_genus_point_one_filtering_df <- GR_by_genus_point_one_zerofree %>%  #create a new dataframe by filtering the old one
  mutate(count_of_nonzeros = rowSums(.[10]!=0)) %>% #create a new column, count_of_nonzeroes, that counts # nonzero values in mean_GR columns
  group_by(famgen_col) %>% #group by the column of interest
  summarize(n = sum(count_of_nonzeros)) %>% #sums the count_of_nonzeros values within each famgen and saves the value as n
  mutate(n) #add n as a new column in the dataframe

#create a new dataframe that merges the old and new dataframes
GR_by_genus_point_one_filter_added <- dplyr::inner_join(GR_by_genus_point_one_zerofree, GR_by_genus_point_one_filtering_df, by= "famgen_col")

#now, remove all rows where n (the name of the count column) is less than 40 (50% of the total number of samples)
GR_by_genus_point_one_filtered <- filter(GR_by_genus_point_one_filter_added, n >= 2)
GR_by_genus_filtered_unique <- unique(GR_by_genus_point_one_filtered$famgen_col) #count unique taxa; as of 05.02.22, we have 54.
GR_by_genus_filtered_unique

#rm(GR_by_genus_point_one, GR_by_genus_point_one_filter_added, GR_by_genus_point_one_zerofree, GR_by_genus_point_one_filtering_df)

#########Step 4: filter to only include paired taxa (i.e. taxa which appear in both tumor and normal samples)
#create a new dataframe where only taxa with growth rates in BOTH tumor and normal samples (i.e. Patient_Blind_ID appears twice) are included
#genus_point_one_paired <- GR_by_genus_point_one_filtered   %>% unite("sorting_col", famgen_col:Patient_Blind_ID, remove = FALSE)
#genus_point_one_paired_only <- subset(genus_point_one_paired,duplicated(sorting_col) | duplicated(sorting_col, fromLast=TRUE))
#after all this filtering, we get 1568 taxa-Patient_Blind_ID combinations. Not sure how this compares to previous results
#I'm not sure why... we're doing this. By definition the data are already paired. 
#genus_point_one_paired_only <- copy(GR_by_genus_point_one_filtered)
#oh, I think there was a mix-up between wide and longform data here.

##########Step 5: split into multiple dataframes by genus and do stats on each genus
#first, edit the dataframe
#drop a couple of columns
#genus_point_one_paired_only <- select(GR_by_genus_point_one_filtered, -c(mean_GR_difference, n, Description))
#rename a couple of columns to make life easier later. Yeah, yeah, I know, this is the reverse of what I did before.
#genus_point_one_paired_only <- genus_point_one_paired_only %>% rename(normal= mean_genus_GR_normal, tumor= mean_genus_GR_tumor) #rename these columns

#colnames(genus_point_one_paired_only)

#convert to long form
#genus_point_one_paired_only <- reshape2::melt(data= genus_point_one_paired_only,
                                    #id.vars= c("Kingdom", "Phylum", "Class", "Order", "famgen_col", "Patient_Blind_ID", "Family", "Genus",
                                    #"Age", "Diagnosis", "TNM", "Statium", "Site"),
                                    #variable.name = "Description",
                                    #value.name = "mean_genus_GR")
#split up by different genera and see where we get. 
split_0.1_GR_genus <- split(GR_by_genus_point_one_filtered, with(GR_by_genus_point_one_filtered, interaction(famgen_col)), drop = TRUE)
#let's see if any taxa have differential GRs among patient samples
genus_0.1_stats <- lapply(split_0.1_GR_genus, function(df){kruskal.test(mean_genus_GR~Isolate, data=df)})

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
#rm(elem, genus_0.1_stats, genus_point_one_paired_only, split_0.1_GR_genus, 
   #genus_0.1_pvals, genus_0.1_qvals, genus_ids_point_one, new_value)

###########Step 8: pull out significant p-values only
genus_0.1_statistics <- genus_0.1_statistics[order(genus_0.1_statistics$genus_0.1_pvals),] #re-order by p-value
#top_ten_genera <- genus_0.1_statistics %>% top_n(-20, genus_0.1_pvals) #note the -20 here to select the smallest p-values (20 taxa)

#no good p-values. 

##########Step 9: save results as a csv file
#write.csv(as.data.frame(genus_0.1_statistics), file="11.08.21_filtered_microbial_GR_all_results_Niccolai2020_data.csv")
#new: filter by sig p-values and only keep significantly different p-values
#genus_0.1_statistics_sig <- genus_0.1_statistics[genus_0.1_statistics$genus_0.1_pvals <= 0.05, ]
#genus_0.1_statistics_qsig <- genus_0.1_statistics[genus_0.1_statistics$genus_0.1_qvals <= 0.05, ]
#write.csv(as.data.frame(genus_0.1_statistics_sig), file="11.08.21.21_filtered_microbial_GR_SIG_results_Niccolai2020_data.csv")

#I guess just save the output.



#DID NOT DO GRAPHING HERE
##########Step 10: GRAPH??!!!
#GR_by_genus_wide <- GR_by_genus_point_one_filtered %>%
  mutate(Difference_direction = case_when(
    mean_GR_difference > 0 ~ "tumor > normal",
    mean_GR_difference < 0 ~ "tumor < normal",
    mean_GR_difference == 0 ~ "Zero change",
  ))

#GRAPH!!!

#count number of samples with Fuso:
freqs <- table(GR_by_genus_wide$Genus)
freqs

#Fusobacterium
ggpaired(subset(GR_by_genus_wide, Genus %in% c("Fusobacterium")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Fusobacterium growth rates, Niccolai2020 data"), subtitle = expression("paired Wilcoxan, p=0.00000230, q= 0.0000827")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.6)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Ruminococcus2
ggpaired(subset(GR_by_genus_wide, Genus %in% c("Ruminococcus2")), cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Ruminococcus growth rates, Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.00182, q= 0.0327")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.4)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

