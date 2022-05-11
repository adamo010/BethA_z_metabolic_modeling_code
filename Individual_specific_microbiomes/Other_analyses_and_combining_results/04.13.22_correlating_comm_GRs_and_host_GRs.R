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

#the goal is to correlate microbiome community and host growth rates. I don't know how this is going to go. 

#data are generated here: 02.22.22_Burns2022_CORDA_model_analyzing_EUstd_GRs.R (host GRs)
#and 09.09.21_Burns2015_comm_GR_combining_files_and_adding_metadata.py (micom GRs)

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Other_analyses_and_combining_results/") 

##########Import data and adjust metadata
#import micom grs
micom_grs <- read.csv(file="/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/09.09.21_Burns2015_V30MSIrun_community_GRs_combined_plus_ranks.csv")
micom_grs <- select(micom_grs, -c(X))
micom_grs <- micom_grs %>% rename(micom_sample_id = Sample_ID, micom_gr = comm_gr) #rename columns for merging
#import host grs
host_grs <- read.csv(file="/Users/adamo010/Documents/CORDA_CRC_human_models/EUstd_diet_model_outputs/02.22.22_Burns2015_human_models_EUstd_GRs_combined.csv", header=FALSE)
host_grs <-host_grs %>% rename(host_sample_id = V1, host_gr = V2)
#import and clean up metadata
metadata= read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/CRC_metadata_rnaseq_edited_by_BA.csv")
metadata_slimmed <- select(metadata, c(Tissue_Tube_ID, Patient_Blind_ID, SampleID, Description)) #drop the columns called ...1 (the index column, and the sample ID column)
rm(metadata)
metadata_slimmed <-metadata_slimmed %>% rename(host_sample_id = Tissue_Tube_ID, micom_sample_id = SampleID)
#merge metadata with host grs
host_grs_plus_metadata <-left_join(host_grs, metadata_slimmed, by = c("host_sample_id" = "host_sample_id"))
#merge in micom data
all_grs_plus_metadata <- left_join(host_grs_plus_metadata, micom_grs, by = c("micom_sample_id"="micom_sample_id"))
#set factors and clean up 
all_grs_plus_metadata <- select(all_grs_plus_metadata, -c(Patient_Blind_ID.x, Description.x, rank, Site, MSI_status, Stage))
all_grs_plus_metadata <- all_grs_plus_metadata %>% rename(Patient_Blind_ID = Patient_Blind_ID.y, Description = Description.y)

#############Correlate!
#setting CI to true; makes things take (a second or two) longer to run, but CI output is nice.
host_micom_gr_corrs <- cor.test(all_grs_plus_metadata$host_gr, all_grs_plus_metadata$micom_gr, use= "pairwise", method = "spearman")
host_micom_gr_corrs #haha fucking nothing.

#############Graph
ggplot(all_grs_plus_metadata, aes(x=host_gr, y=micom_gr)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("Host-microbiome growth rate correlations"), subtitle= expression("Spearman correlation p = 0.7055")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="Microbiome growth rate, mmmol/(gDW * h)") +
  scale_y_continuous(name = "Host growth rate, mmmol/(gDW * h)") 


