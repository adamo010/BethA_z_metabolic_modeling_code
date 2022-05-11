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
library("readxl")
library("forcats")
library("viridis")

#the goal here is to analyze the MICOM indiv spp growth rates from communities grown on host spent medium
#start with 11.08.21_Niccolai2020_indiv_spp_GR_analysis and 03.21.22_combining_host_and_microbiomes_EUstd_host_GR_analysis.R

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_analyses/") 

##########step 1: import data and clean up
#this dataframe contains the longform data; see 08.31.21_Hale2018_comm_GR_combining_files_and_adding_metadata.py for generating it. 
indiv_GRs_all <- read_csv("04.04.22_combining_host_and_microbiomes_EUstd_microbiome_indiv_spp_GRs_collapsed_by_genus.csv")
indiv_GRs_all <- select(indiv_GRs_all, -c(...1)) #drop the columns called ...1 (the index column, and the sample ID column)
#set theme to have nice colours
theme_set(theme_bw())
#set Patient_Blind_ID as factor
indiv_GRs_all$Patient_Blind_ID <- as.factor(indiv_GRs_all$Patient_Blind_ID)

#########Step 2: reorder and add columns for filtering and stats
#add a new column, famgen, for sorting purposes- some genera have the same name but different families. 
GR_by_genus_point_one <- indiv_GRs_all  %>% unite("famgen_col", family:genus, remove = FALSE) 
#drop useless columns
#GR_by_genus_point_one = subset(GR_by_genus_point_one, select = -c(...1, index, host_subject_id) )
colnames(GR_by_genus_point_one) #check
#now rearrange so famgen and Patient_Blind_ID are next to each other. 
GR_by_genus_point_one <- GR_by_genus_point_one[, c(1,2,3,4,5,6,10,7,8,9,11,12,13,14,15,16)] #this is edited from previous version, as my output is a bit different
#how many unique taxa are we including here? 
GR_by_genus_unique <- unique(GR_by_genus_point_one$famgen_col) #84! ok

#########Step 3 NEW: filter by taxa which appear in at least 20% of the total number of samples 
#first, flip our dataframe from long to wide. Need mean_GR_differences sorted by combo_type
#drop some columns to allow proper spread
GR_by_genus_point_one_wide <- select(GR_by_genus_point_one, -c("host_model_description", "host_model_id", "microbiome_model_description",
                                                               "microbiome_model_id", "host_plus_microbiome_model"))
GR_by_genus_point_one_wide <- spread(GR_by_genus_point_one_wide, pair_type, mean_genus_GR)
GR_by_genus_point_one_wide[is.na(GR_by_genus_point_one_wide)] = 0 #replace all nans with 0s
colnames(GR_by_genus_point_one_wide) #check

GR_by_genus_point_one_wide <- GR_by_genus_point_one_wide %>%# Compute row sums
  mutate(sum = rowSums(GR_by_genus_point_one_wide[, c(10:13)] != 0))
#now we have a new column, sum, that counts the number of non-zero GRs
#now, remove all rows (OTUs) where all GRs are 0
GR_by_genus_point_one_zerofree <- GR_by_genus_point_one_wide[(GR_by_genus_point_one_wide$sum !=0),]

#how many unique taxa are we including here? 55. Not bad. 
GR_by_genus_zerofree_unique <- unique(GR_by_genus_point_one_zerofree$famgen_col)

colnames(GR_by_genus_point_one_zerofree) #check column names; will be used below
#what do we actually want here? want to limit the number of taxa that we have to look at. 
#first, we need at least two groups to compare. So, need GR_by_genus_point_one_mintwo, where sum >=2
GR_by_genus_point_one_mintwo <- GR_by_genus_point_one_wide[(GR_by_genus_point_one_wide$sum >=2),]
GR_by_genus_mintwo_unique <- unique(GR_by_genus_point_one_mintwo$famgen_col)
#that gets us 43 taxa.
#now, we need each taxon to appear in... at least 3 PBIDs? That's enough for stats, right?
#uh, maybe not for multiple test corrections? maybe need seven samples, b/c we're doing six tests?
#let's start with three PBIDs, and see where we get.

# create a filtering dataframe that contains counts of each PBID within famgen 
#can we just count number of rows for each famgen_col? Since we converted to wide form, each PBID-taxon is a separate, unique row
GR_by_genus_point_one_filtered <- GR_by_genus_point_one_mintwo %>%  #create a new dataframe by filtering the old one
  group_by(famgen_col) %>% #group by the column of interest
  mutate(count=n()) #add n as a new column in the dataframe

#create a new dataframe that merges the old and new dataframes
#GR_by_genus_point_one_filter_added <- dplyr::inner_join(GR_by_genus_point_one_zerofree, GR_by_genus_point_one_filtering_df, by= "famgen_col")

#now, remove all rows where n (the name of the count column) is less than 50% of the total number of samples.
#NOTE here: max # counts should be 88*4 = 352. B/c there are four possible combos: TN, NT, NN, TT
#so, half max is 176. I'm not sure if this is exactly the right approach, but we'll see what it gets us.
#HAHAHA 176 got us ONE sample. What's fuso at? very low. 

#here, count is the number of occurances of that famgen. Let's go with a cutoff of >=7
GR_by_genus_point_one_filtered_highfreq <- filter(GR_by_genus_point_one_filtered, count >= 7)
GR_by_genus_filtered_unique <- unique(GR_by_genus_point_one_filtered_highfreq$famgen_col) #count unique taxa
GR_by_genus_filtered_unique #now there are 29. 
rm(GR_by_genus_point_one, GR_by_genus_point_one_filter_added, GR_by_genus_point_one_zerofree, GR_by_genus_point_one_filtering_df)
#step coming up anyway. 

#########Step 4: filter to only include paired taxa (i.e. taxa which appear in both tumor and normal samples)
#let's just... okay, the only way to do stats is to have at least two GR_by_genus_point_one_zerofree. Which we do, now.
#GR_by_genus_point_one_filtered_highfreq has only OTUs which appear at least twice within a PBID, and only OTUs which appear
#in at least 7 PBIDs. & was chosen to allow good stats.

##########Step 5: split into multiple dataframes by genus and do stats on each genus
#first, edit the dataframe
#drop a couple of columns
genus_point_one_paired_only <- select(GR_by_genus_point_one_filtered_highfreq, -c(sum, count, OTU_ID))
colnames(genus_point_one_paired_only)

#convert to long form
genus_point_one_paired_only <- reshape2::melt(data= genus_point_one_paired_only,
                                              id.vars= c("kingdom", "phylum", "class", "order", "famgen_col", "family", "genus", "Patient_Blind_ID"),
                                              variable.name = "Description",
                                              value.name = "mean_genus_GR")
#split up by different genera and see where we get. 29 genera
split_0.1_GR_genus <- split(genus_point_one_paired_only, with(genus_point_one_paired_only, interaction(famgen_col)), drop = TRUE)

#the below code takes all the dataframes (one per taxon) in split_0.1_GR_genus and applies an ANOVA test.
genus_0.1_stats <- lapply(split_0.1_GR_genus, function(df){aov(mean_genus_GR~Description, data=df)})
#from the ANOVA test, we pull out summary (gives P-values) and Tukey HSD tests
genus_0.1_summary <- lapply(genus_0.1_stats, function(df){summary(df)})
genus_0.1_hsd <- lapply(genus_0.1_stats, function(df){TukeyHSD(df)})


##########Step 6: statistics
#make lists of p values
genus_0.1_pvals <- c() #initialize list
for (elem in genus_0.1_summary) {
  new_value = elem[[1]][["Pr(>F)"]] #extract p-value from the dataframe
  new_value2 = new_value[1] #trim off weird nans
  genus_0.1_pvals <-c(genus_0.1_pvals, new_value2) #add p-value to genus_0.1_pvals
}
rm(new_value, new_value2)

#make a list of q-values
genus_0.1_qvals <- p.adjust(genus_0.1_pvals, method = "fdr")

#make a list of genus_ids
genus_ids_point_one <- c()
for(elem in split_0.1_GR_genus){
  new_value = elem$famgen_col[1]
  genus_ids_point_one <- c(genus_ids_point_one, new_value)}
#merge all lists together. 
genus_0.1_statistics <- data.frame(genus_ids_point_one, genus_0.1_pvals, genus_0.1_qvals)

genus_0.1_hsds <- c() #initialize list
for (elem in genus_0.1_hsd) {
  new_value = elem[[1]][["Pr(>F)"]] #extract p-value from the dataframe
  new_value2 = new_value[1] #trim off weird nans
  genus_0.1_hsds <-c(genus_0.1_pvals, new_value2) #add p-value to genus_0.1_pvals
}
rm(new_value, new_value2)

#okay, now we extract useful info from hsd results
hsd_res_df <- data.frame(matrix(ncol = 6)) #create empty dataframe
colnames(hsd_res_df) <- c('combo_type', 'diff', 'lwr', 'upr', 'p adj', 'famgen') #name columns in dataframe
elem_counter = 0 #initialize counter
for (elem in genus_0.1_hsd) {
  elem_counter = elem_counter+1 #add counter
  df <- elem$Description #extract useful info from each element in the genus_0.1_hsd list of results
  df2 <- data.table::as.data.table(df, keep.rownames=TRUE) #convert output to a data table, keeping row names
  df2 <- rename(df2, combo_type=rn) #rename column
  new_value = names(genus_0.1_hsd[elem_counter]) #get combo_type name
  df2 <- df2 %>% mutate(famgen = new_value) #add combo_type name as constant value in famgen column
  hsd_res_df <- rbind(hsd_res_df, df2) #append this element's results to hsd_res_df 
}
rm(df,df2,elem,elem2,elem_counter,new_value)


##########Step 9: save results as a csv file. START HERE.
write.csv(as.data.frame(genus_0.1_statistics), file="04.07.22_combining_host_and_microbiomes_EUstd_microbiome_indiv_spp_GRs_ANOVA_results_statistics.csv")
write.csv(as.data.frame(hsd_res_df), file="04.07.22_combining_host_and_microbiomes_EUstd_microbiome_indiv_spp_GRs_Tukey_hsd_results_statistics.csv")
write.csv(as.data.frame(genus_point_one_paired_only), file="04.07.22_combining_host_and_microbiomes_EUstd_microbiome_indiv_spp_GRs_longform_data.csv")

#########Step 10: graph. 
y_comparisons <- list( c("tumor_tumor", "tumor_normal"), c("tumor_tumor", "normal_tumor"), c("tumor_tumor", "normal_normal"),
                       c("tumor_normal", "normal_tumor"), c("tumor_normal", "normal_normal"), c("normal_tumor", "normal_normal"))
#graph: Blautia
ggplot(subset(genus_point_one_paired_only, genus %in% c("Blautia")), aes(x= Description, y= mean_genus_GR, color=Description)) +
  geom_boxplot(aes(x= Description, y= mean_genus_GR)) + 
  geom_dotplot(binaxis='y', stackdir='center', stackratio=0.75, dotsize=0.6, aes(x= Description, y= mean_genus_GR, fill=Description)) +
  labs(title = expression("MICOM Blautia growth rates, Burns2015 data"), 
       subtitle= expression("p-values from uncorrected pairwise Mann-Whitney U tests")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), aspect.ratio=1.5/2) +
  scale_y_continuous(name = "Community growth\nrate, mmol/(gDW * h)", limits=c(-0.1,0.6)) +
  scale_x_discrete(labels = c('Normal\nmicrobiome,\nnormal host','Normal\nmicrobiome,\ntumor host',
                              'Tumor\nmicrobiome,\nnormal host','Tumor\nmicrobiome,\ntumor host'),
                   name= "Patient microbiome and host tissue descriptions") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "none") +
  guides(fill='none') +
  stat_compare_means(comparisons= y_comparisons, method = "wilcox.test", size= 3, label.y = c(0.40,0.44,0.48,0.52,0.56,0.60))

#another way of looking at the data

ggplot(subset(genus_point_one_paired_only, genus %in% c("Blautia")), aes(x= Patient_Blind_ID, y= mean_genus_GR, color=Description, fill= Description)) +
  geom_dotplot(binaxis='y', stackdir='center', stackratio=0.75, dotsize=0.6, aes(x= Patient_Blind_ID, y= mean_genus_GR))+
  labs(title = expression("MICOM Blautia growth rates by\nPatient_Blind_ID, Burns2015 data"), color= "Host-microbiome combination", 
       fill= "Host-microbiome combination") +
  scale_y_continuous(name = "Microbial growth\nrate, mmol/(gDW * h)", limits=c(0,0.4)) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), aspect.ratio=3.5/4) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2")


