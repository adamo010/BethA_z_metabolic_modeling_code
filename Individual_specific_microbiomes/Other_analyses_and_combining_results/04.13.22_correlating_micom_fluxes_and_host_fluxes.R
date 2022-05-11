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

#the goal is to correlate microbiome and host exchange fluxes. 

#data are generated here: 02.22.22_Burns2022_CORDA_model_analyzing_EUstd_GRs.R (host GRs)
#and 09.09.21_Burns2015_comm_GR_combining_files_and_adding_metadata.py (micom GRs)

#analyses borrowed from 4.12.22_correlating_Burns_GEs_and_all_microbiome_metrics.R and 02.17.22_Burns2022_CORDA_model_analyzing_EUstd_fluxes_exchanges_only.R

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Other_analyses_and_combining_results/") 

############step 1: import files and adjust factors
micom_flux_values <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/10.25.21_filtered_sum_exchange_fluxes_only_alldata_Burns2015_data.csv") #used to be fluxes_top_diffs
micom_flux_values <- select(micom_flux_values, -c(...1))
micom_flux_stats <-read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/10.25.21_filtered_sum_exchange_fluxes_only_stats_Burns2015_data.csv") #used to be fluxes_all
micom_flux_stats <- select(micom_flux_stats, -c(...1))

host_flux_values <- read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/EUstd_diet_model_outputs/02.17.22_CORDA_flux_exchanges_alldata_Burns2015_data_EUstd.csv")
host_flux_values <- select(host_flux_values, -c(...1))
host_flux_stats <- read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/EUstd_diet_model_outputs/02.17.22_CORDA_flux_exchanges_stats_Burns2015_data_EUstd.csv")
host_flux_stats <- select(host_flux_stats, -c(...1))

############step 2: filter by top 20 significantly different exchange pathways. .
#note: here we're filtering by p-values, but the order is the same as if we filtered by p-value
micom_flux_stats <- micom_flux_stats %>% top_n(-20, by_flux_pvals)
host_flux_stats <- host_flux_stats %>% top_n(-20, by_flux_pvals)
#take a vector of the fluxpath_ids
top_micom_flux_names <- pull(micom_flux_stats, flux_ids) #make a vector from a column
top_host_flux_names <- pull(host_flux_stats, fluxpath_name)
#use that to filter on the new column
top_micom_flux_data <- filter(micom_flux_values, fluxpath_name %in% top_micom_flux_names)
top_host_flux_data <- filter(host_flux_values, fluxpath_name %in% top_host_flux_names)
#then, rename columns
top_micom_flux_data <- top_micom_flux_data %>% rename(fluxes_normal = normal, fluxes_tumor = tumor) #rename columns
top_micom_flux_data$Patient_Blind_ID <- as.factor(top_micom_flux_data$Patient_Blind_ID)
top_host_flux_data <- top_host_flux_data %>% rename(fluxes_normal = normal, fluxes_tumor = tumor) #rename columns
top_host_flux_data$Patient_Blind_ID <- as.factor(top_host_flux_data$Patient_Blind_ID)
#clean up
rm(host_flux_values, micom_flux_values, top_host_flux_names, top_micom_flux_names)

##########step 3: calculate logfold changes... I think I'm going to pass on this for now. Taking the log of flux data can be a mess.

############Step 4: drop cols, convert dataframes to wide, and set row names as Patient_Blind_ID for correlating
#Do both tumor and normal samples here
top_micom_flux_data_wide_normal <- subset(top_micom_flux_data, select = -c(fluxes_tumor,diff,Difference_direction))
top_micom_flux_data_wide_normal <- spread(top_micom_flux_data_wide_normal, fluxpath_name, fluxes_normal) #convert dataframe to wide
top_micom_flux_data_wide_normal <- column_to_rownames(top_micom_flux_data_wide_normal, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names
top_micom_flux_data_wide_tumor <- subset(top_micom_flux_data, select = -c(fluxes_normal,diff,Difference_direction))
top_micom_flux_data_wide_tumor <- spread(top_micom_flux_data_wide_tumor, fluxpath_name, fluxes_tumor) #convert dataframe to wide
top_micom_flux_data_wide_tumor <- column_to_rownames(top_micom_flux_data_wide_tumor, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names

top_host_flux_data <- select(top_host_flux_data, Patient_Blind_ID, fluxpath_name, fluxes_normal, fluxes_tumor) #sort dataframe columns. 
top_host_flux_data <- top_host_flux_data[order(top_host_flux_data$Patient_Blind_ID),]
top_host_flux_data_wide_normal <- subset(top_host_flux_data, select = -c(fluxes_tumor))
top_host_flux_data_wide_normal <- spread(top_host_flux_data_wide_normal, fluxpath_name, fluxes_normal) #convert dataframe to wide
top_host_flux_data_wide_normal <- column_to_rownames(top_host_flux_data_wide_normal, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names
top_host_flux_data_wide_tumor <- subset(top_host_flux_data, select = -c(fluxes_normal))
top_host_flux_data_wide_tumor <- spread(top_host_flux_data_wide_tumor, fluxpath_name, fluxes_tumor) #convert dataframe to wide
top_host_flux_data_wide_tumor <- column_to_rownames(top_host_flux_data_wide_tumor, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names

##############CORRELATING##############
##########Step 1: Spearman correlations with Benjamini-Hochberg corrections
#setting CI to true; makes things take (a second or two) longer to run, but CI output is nice.
normal_micom_host_flux_correlations <- corr.test(top_host_flux_data_wide_normal, top_micom_flux_data_wide_normal,
                                                 use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5)
tumor_micom_host_flux_correlations <- corr.test(top_host_flux_data_wide_tumor, top_micom_flux_data_wide_tumor,
                                                 use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5)
###########Step 2: get results into dataframe
#corrected and uncorrected p-values, and r-values
normal_correlations_corrected <- as.data.frame(normal_micom_host_flux_correlations$p.adj) #this makes a dataframe
normal_correlations_corrected <- tibble::rownames_to_column(normal_correlations_corrected, "host_fluxpath") #THANKS DPLYR, hero of R!
normal_correlations_uncorrected <- as.data.frame(normal_micom_host_flux_correlations$p) #this makes a dataframe
normal_correlations_uncorrected <- tibble::rownames_to_column(normal_correlations_uncorrected, "host_fluxpath") #THANKS DPLYR, hero of R!
normal_correlations_rvals <- as.data.frame(normal_micom_host_flux_correlations$r) #this makes a dataframe
normal_correlations_rvals <- tibble::rownames_to_column(normal_correlations_rvals, "host_fluxpath") #THANKS DPLYR, hero of R!

tumor_correlations_corrected <- as.data.frame(tumor_micom_host_flux_correlations$p.adj) #this makes a dataframe
tumor_correlations_corrected <- tibble::rownames_to_column(tumor_correlations_corrected, "host_fluxpath") #THANKS DPLYR, hero of R!
tumor_correlations_uncorrected <- as.data.frame(tumor_micom_host_flux_correlations$p) #this makes a dataframe
tumor_correlations_uncorrected <- tibble::rownames_to_column(tumor_correlations_uncorrected, "host_fluxpath") #THANKS DPLYR, hero of R!
tumor_correlations_rvals <- as.data.frame(tumor_micom_host_flux_correlations$r) #this makes a dataframe
tumor_correlations_rvals <- tibble::rownames_to_column(tumor_correlations_rvals, "host_fluxpath") #THANKS DPLYR, hero of R!

###########Step 3: get data into longform
normal_correlations_corrected_long <- normal_correlations_corrected %>%
  melt(id.var="host_fluxpath") %>%
  arrange(host_fluxpath, variable) 
normal_correlations_corrected_long <- rename(normal_correlations_corrected_long, c("microbiome_fluxpath"="variable", "p_value_adj"="value")) #rename columns
normal_correlations_uncorrected_long <- normal_correlations_uncorrected %>%
  melt(id.var="host_fluxpath") %>%
  arrange(host_fluxpath, variable) 
normal_correlations_uncorrected_long <- rename(normal_correlations_uncorrected_long, c("microbiome_fluxpath"="variable", "p_value_raw"="value")) #rename columns
normal_correlations_rvals_long <- normal_correlations_rvals %>%
  melt(id.var="host_fluxpath") %>%
  arrange(host_fluxpath, variable) 
normal_correlations_rvals_long <- rename(normal_correlations_rvals_long, c("microbiome_fluxpath"="variable", "R_value"="value")) #rename columns

tumor_correlations_corrected_long <- tumor_correlations_corrected %>%
  melt(id.var="host_fluxpath") %>%
  arrange(host_fluxpath, variable) 
tumor_correlations_corrected_long <- rename(tumor_correlations_corrected_long, c("microbiome_fluxpath"="variable", "p_value_adj"="value")) #rename columns
tumor_correlations_uncorrected_long <- tumor_correlations_uncorrected %>%
  melt(id.var="host_fluxpath") %>%
  arrange(host_fluxpath, variable) 
tumor_correlations_uncorrected_long <- rename(tumor_correlations_uncorrected_long, c("microbiome_fluxpath"="variable", "p_value_raw"="value")) #rename columns
tumor_correlations_rvals_long <- tumor_correlations_rvals %>%
  melt(id.var="host_fluxpath") %>%
  arrange(host_fluxpath, variable) 
tumor_correlations_rvals_long <- rename(tumor_correlations_rvals_long, c("microbiome_fluxpath"="variable", "R_value"="value")) #rename columns

##########Step 4: merge and save dataframes
normal_correlations_long_merged <- merge(merge(normal_correlations_corrected_long,
                                               normal_correlations_uncorrected_long, all=TRUE),
                                         normal_correlations_rvals_long, all=TRUE)
normal_correlations_long_merged <- normal_correlations_long_merged[order(normal_correlations_long_merged$p_value_raw),] #sort by p-value
tumor_correlations_long_merged <- merge(merge(tumor_correlations_corrected_long,
                                              tumor_correlations_uncorrected_long, all=TRUE),
                                        tumor_correlations_rvals_long, all=TRUE)
tumor_correlations_long_merged <- tumor_correlations_long_merged[order(tumor_correlations_long_merged$p_value_raw),] #sort by p-value

#pull in host metadata
host_flux_metadata <- read.table(file="/Users/adamo010/Documents/CORDA_CRC_human_models/recon_model_fluxpaths.tsv", sep = '\t', header = TRUE)
host_flux_metadata <- select(host_flux_metadata, c(abbreviation, description))
host_flux_metadata <- rename(host_flux_metadata, c("host_fluxpath_description"="description"))
host_flux_metadata_exchanges <- dplyr::filter(host_flux_metadata, grepl("EX_",abbreviation)) #extract exchange reactions
host_flux_metadata_exchanges <- host_flux_metadata_exchanges %>%
  mutate_at("abbreviation", str_replace, "EX_", "")
host_flux_metadata_exchanges$abbreviation <- str_sub(host_flux_metadata_exchanges$abbreviation, 1, str_length(host_flux_metadata_exchanges$abbreviation)-3)

#re-merge
normal_correlations_long_merged <- left_join(normal_correlations_long_merged, host_flux_metadata_exchanges, by=c("host_fluxpath"="abbreviation"))
tumor_correlations_long_merged <- left_join(tumor_correlations_long_merged, host_flux_metadata_exchanges, by=c("host_fluxpath"="abbreviation"))

#pull in microbiome metadata
micom_flux_metadata = read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Common_growth_medium_runs/May_June2021_MICOM_flux_and_pathway_analyses/Microbiome_flux_key_full_June2021.csv")
micom_flux_metadata <- select(micom_flux_metadata, c(abbreviation, description))
micom_flux_metadata <- rename(micom_flux_metadata, c("microbiome_fluxpath_description"="description"))
micom_flux_metadata_exchanges <- dplyr::filter(micom_flux_metadata, grepl("EX_",abbreviation)) #extract exchange reactions

#re-merge
normal_correlations_long_merged <- left_join(normal_correlations_long_merged, micom_flux_metadata, by=c("microbiome_fluxpath"="abbreviation"))
tumor_correlations_long_merged <- left_join(tumor_correlations_long_merged, micom_flux_metadata_exchanges, by=c("microbiome_fluxpath"="abbreviation"))

#########Step 5: save the files, just in case
write.csv(as.data.frame(normal_correlations_long_merged), file="04.13.22_Burns2015_micom_by_host_flux_correlations_normal_samples.csv")
write.csv(as.data.frame(tumor_correlations_long_merged), file="04.13.22_Burns2015_micom_by_host_flux_correlations_tumor_samples.csv")

#clean up
rm(host_flux_stats, micom_flux_stats, normal_correlations_corrected, normal_correlations_corrected_long, normal_correlations_rvals,normal_correlations_rvals_long)
rm(normal_correlations_uncorrected, normal_correlations_uncorrected_long, normal_micom_host_flux_correlations, top_host_flux_data, top_micom_flux_data)
rm(tumor_correlations_corrected,tumor_correlations_corrected_long, tumor_correlations_rvals, tumor_correlations_rvals_long, tumor_correlations_uncorrected)
rm(tumor_correlations_uncorrected_long,tumor_micom_host_flux_correlations)

#########Step 6: prepare for graphing
#filter to only keep rows where unadjusted P<0.05
normal_correlations_long_merged_sigonly <- filter(normal_correlations_long_merged, p_value_raw < 0.05) #0.05 gets us 37 correlations
tumor_correlations_long_merged_sigonly <- filter(tumor_correlations_long_merged, p_value_raw < 0.05) #0.05 gets us 8 correlations

#save the files, just in case
write.csv(as.data.frame(normal_correlations_long_merged_sigonly), file="04.13.22_Burns2015_micom_by_host_flux_correlations_normal_samples_sigto0.01_correlations.csv")
write.csv(as.data.frame(tumor_correlations_long_merged_sigonly), file="04.13.22_Burns2015_micom_by_host_flux_correlations_tumor_samples_sigto0.01_correlations.csv")

#collect data to graph.
normal_fluxpaths_combined <-cbind(top_host_flux_data_wide_normal, top_micom_flux_data_wide_normal) #NOTE: we can only do the cbind b/c we have the same rows in the same order.
tumor_fluxpaths_combined <-cbind(top_host_flux_data_wide_tumor, top_micom_flux_data_wide_tumor) #NOTE: we can only do the cbind b/c we have the same rows in the same order.

#########Step 7: graphing
#for graphing, data comes from fluxpaths_and_ges_combined, and p-values come from fluxpath_correlations_all_uncorrected_long
#first, need to remove (e) from all column headers
colnames(normal_fluxpaths_combined)<-gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(normal_fluxpaths_combined)))
colnames(tumor_fluxpaths_combined)<-gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(tumor_fluxpaths_combined)))
#I am aware that this doesn't specifically remove (e), but it's the only way to get rid of the ellipses AND it gets rid of e, so...

#X axis is microbiome, y axis is host
ggplot(normal_fluxpaths_combined, aes(x=EX_thr_L, y=duri)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("Host deoxyiridine exchange flux vs. microbiome\nL-Threonine exchange flux, normal samples")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.margin=unit(c(7,1,1,1),"mm")) +
  guides(fill="none") +
  scale_x_continuous(name ="Microbiome L-Threonine exchange flux, mmmol/(gDW * h)") +
  scale_y_continuous(name = "Host deoxyuridine exchange flux, mmmol/(gDW * h)") +
  annotate(geom="text", x=-4000, y=800, label="Spearman correlation\np.adj = 0.0787",color="black")

ggplot(tumor_fluxpaths_combined, aes(x=EX_lac_D, y=	leuktrE4)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("Host Leukotriene E4 exchange flux vs. microbiome\nD-lactate exchange flux, tumor samples")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.margin=unit(c(7,1,1,1),"mm")) +
  guides(fill="none") +
  scale_x_continuous(name ="Microbiome D-lactate exchange flux, mmmol/(gDW * h)") +
  scale_y_continuous(name = "Host Leukotriene E4 exchange flux, mmmol/(gDW * h)") +
  annotate(geom="text", x=2000, y=-750, label="Spearman correlation\np.raw = 0.00780",color="black")





