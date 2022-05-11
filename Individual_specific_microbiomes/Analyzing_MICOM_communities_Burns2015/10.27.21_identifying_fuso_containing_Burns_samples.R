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

#the goal here is to figure out which SampleIDs contain Fusobacterium. i need these to compare gene expression between +/- fuso samples
#this seemed like as good a place as any to start. 

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/") 

##########Step 1: import data and adjust factors
#It is important to note here that I have data (and have done analyses) for both family and genus level. However, I went with genus level
#moving forward because... most of the family level associations are basically due to a single genus anyway, I think is what I decided. 
#NOTE: this is a temporary file; edit as needed
GR_by_genus_point_one <- read_csv("09.09.21_Burns2015_indiv_spp_GRs_collapsed_to_genus.csv")
GR_by_genus_point_one$Patient_Blind_ID <- as.factor(GR_by_genus_point_one$Patient_Blind_ID)
#colnames(GR_by_genus_point_one) #print column names as needed

########Step 2: create a fuso-only Df and filter by tumor samples with growing fuso
fuso_only_samples <- filter(GR_by_genus_point_one, genus == "Fusobacterium")
fuso_tumor_samples <- filter(fuso_only_samples, Description == "tumor")
fuso_growing_tumor_samples <- filter(fuso_tumor_samples, mean_genus_GR > 0)

########Step 3: get a list of samples with fuso growing
growing_fuso_sampleids <- pull(fuso_growing_tumor_samples, SampleID)

########Step 4: create a new dataframe with all (tumor) sample_ids and a second column, fuso_present, with a value of 'yes' if that sampleid in
#growing_fuso_sampleids, and 'no' if it is not.
new_metadata_df <- select(GR_by_genus_point_one, c(SampleID, Description, Patient_Blind_ID)) #create new metadata df with select columns of interest.
new_metadata_df <- new_metadata_df %>% distinct() #remove duplicates
#I think the easiest way to do this is to split up the df, add my column values, and re-stack the dfs
fuso_pos_metadata_df <- new_metadata_df[new_metadata_df$SampleID %in% growing_fuso_sampleids, ] #create a new df that is only fuso positive samples
fuso_pos_metadata_df <- fuso_pos_metadata_df %>% add_column(Fuso_present = "yes") #add a column called Fuso_present and fill it with 'yes' valyes
fuso_normal_samples_df <- filter(new_metadata_df, Description == "normal")
fuso_normal_samples_df <- fuso_normal_samples_df %>% add_column(Fuso_present = "NA_normal_tissue") 
normal_sampleids <- pull(fuso_normal_samples_df, SampleID)
normal_and_pos_sampleids <- c(normal_sampleids, growing_fuso_sampleids)
fuso_neg_metadata_df <- subset(new_metadata_df, !(SampleID %in% normal_and_pos_sampleids))
fuso_neg_metadata_df <- fuso_neg_metadata_df %>% add_column(Fuso_present = "no")

new_metadata_df2 <- rbind(fuso_pos_metadata_df, fuso_normal_samples_df, fuso_neg_metadata_df)

write.csv(as.data.frame(new_metadata_df2), file="10.28.21_key_for_identifying_Fuso_positive_tumor_samples_Burns2015.csv")

