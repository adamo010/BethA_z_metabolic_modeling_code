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
library("psych")
library("openxlsx")
library("magrittr")
library("purrr")
library("reshape2")

#the goal here is to import differentially expressed genes and significantly different microbial fluxpaths from the Burns RNAseq data 
#and my metabolic models (based on the Burns microbiome data) and correlate them.

#based on 05.28.21_correlating_Burns_125_metabolic_GEs_20_microbial_inoutputs.R
#and Date_TBA_correlating_Burns_125_metabolic_GEs_20_microbial_fluxpaths.R
#and 09.13.21_correlating_Burns_GEs_and_fluxpaths.R
#and 09.17.21_correlating_Burns_GEs_and_flux_subsystems.R 
#and 09.28.21_correlating_Burns_GEs_and_flux_subsystems.R
#and 11.12.21_correlating_Burns_GEs_and_Fuso_fluxes.R

#for the generation of the data used in this script, see
#11.12.21_calculating_log2foldchange_Burns_RNAseq_filtered_add_metabolism_filter_Fuso_samples_only.py (for GEs)
#11.08.21_Burns2015_Fuso_fluxpath_analysis_exchanges_only.R (for fluxes)- NOTE HERE that we're using median values moving forward, instead
#of average values. 

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/") 

##############Fluxpaths##############
top_fluxes <- read_csv("12.20.21_Fusobacterium_exchange_fluxes_only_stats_Burns2015_data.csv")
fluxes_all <-read_csv("12.20.21_Fusobacterium_exchange_fluxes_only_alldata_Burns2015_data.csv")

#check on number of unique PBIDs
unique_PBIDs <- unique(fluxes_all$Patient_Blind_ID)
#all right, 22 here, and 22 from host gene expression

############step 2: filter fluxes_all 
####now, the filtering step. here we decide how many correlations we want to do.
#have two options: use all (43), or use only those with qvals > 0.05 (6) or use only those with pvals < 0.05 (9). Create vectors for each
all_flux_names <- pull(top_fluxes,fluxpath_name) #make a vector from a column:thanks dplyr!
qval_flux_names <- all_flux_names[1:6] #select top 6 elements from the vector
pval_flux_names <- all_flux_names[1:9] #select top 9 elements from the vector

#use that to filter on the new column
fluxes_all_diffs <- filter(fluxes_all, fluxpath_name %in% all_flux_names) #thanks again dplyr!
fluxes_qval_diffs <- filter(fluxes_all, fluxpath_name %in% qval_flux_names) #thanks again dplyr!
fluxes_pval_diffs <- filter(fluxes_all, fluxpath_name %in% pval_flux_names) #thanks again dplyr!

#set Patient_Blind_ID as factor.
fluxes_all_diffs$Patient_Blind_ID <- as.factor(fluxes_all_diffs$Patient_Blind_ID)
fluxes_qval_diffs$Patient_Blind_ID <- as.factor(fluxes_qval_diffs$Patient_Blind_ID)
fluxes_pval_diffs$Patient_Blind_ID <- as.factor(fluxes_pval_diffs$Patient_Blind_ID)

#we're missing a Patient_Blind_ID for qval_diffs
#pffffff, what does this mean..
unique_qids <- unique(fluxes_qval_diffs$Patient_Blind_ID)
missing = setdiff(unique_PBIDs, unique_qids) #find the missing PBID. Remove it from host data later on

############step 3: split dataframes into imports (+ve) and exports (-ve)
#NOT NEEDED HERE

#############step 4: NON ESSENTIAL: take a quick look at the graphs:
#Boxplots!
ggplot(data=fluxes_all_diffs, aes(x=fluxpath_name, y= tumor, fill=fluxpath_name)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  labs(title= "Flux subsystem activity, tumor samples") +
  scale_x_discrete(name ="Fluxpath subsystem name") +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  scale_y_continuous(name = "Flux, mmmol/(gDW * h)") +
  theme(legend.position="none", plot.title = element_text(hjust=0.5, vjust=2, size=12))

#huh. not bad.

############step 5: calculate log2(import) and log2(export) fluxes for tumor samples only- needed b/c we're taking log2fold change of genes.

#calculate log2 values and put into new column
#fluxes_all_diffs[,c("log2_fluxes_tumor")] <- log2(fluxes_all_diffs[,c("tumor")])
#fluxes_all_diffs[,c("log2_fluxes_normal")] <- log2(fluxes_all_diffs[,c("normal")])
#NaNs produced: is this bad? A Problem? Who knows?
#reader, it was a problem. That was not solved by adding a constant 1.
#oh, it's that you can't take the log of a negative value. So, the way to get around this is to....
#standby. need to check original code where I collapsed by subsystems.
#FIXED: couldn't figure out how to account for zeroes (reverse Rx fluxes), so I took medians instead of means within a subsystem.
#kids, this did not work with fuso-only data.
#you know what? No transform. Let's see how taht works.

#fluxes_qval_diffs[,c("log2_fluxes_tumor")] <- log2(fluxes_qval_diffs[,c("tumor")])
#fluxes_qval_diffs[,c("log2_fluxes_normal")] <- log2(fluxes_qval_diffs[,c("normal")])
#
#fluxes_pval_diffs[,c("log2_fluxes_tumor")] <- log2(fluxes_pval_diffs[,c("tumor")])
#fluxes_pval_diffs[,c("log2_fluxes_normal")] <- log2(fluxes_pval_diffs[,c("normal")])

#add a differences column
#fluxes_all_diffs$log2_flux_diffs <- (fluxes_all_diffs$log2_fluxes_tumor - fluxes_all_diffs$log2_fluxes_normal)
#fluxes_qval_diffs$log2_flux_diffs <- (fluxes_qval_diffs$log2_fluxes_tumor - fluxes_qval_diffs$log2_fluxes_normal)
#fluxes_pval_diffs$log2_flux_diffs <- (fluxes_pval_diffs$log2_fluxes_tumor - fluxes_pval_diffs$log2_fluxes_normal)

#another way to do this
#fluxes_top_diffs[,c("log2_flux_diffs_V3")] <- log2(fluxes_top_diffs$fluxes_tumor / fluxes_top_diffs$fluxes_normal)

############Step 6: drop cols, convert dataframes to wide, and set row names as Patient_Blind_ID for correlating
#remove all but three columns
fluxes_all_diffs_cleaned <- subset(fluxes_all_diffs, select = -c(...1, Difference_direction, normal, tumor, fluxpath_description) )
#oh FFS, have to add underscores to flux descriptions or they won't get graphed. 
fluxes_all_diffs_cleaned$fluxpath_name <- gsub(" ", "_", fluxes_all_diffs_cleaned$fluxpath_name)
fluxes_all_diffs_cleaned_wide <- spread(fluxes_all_diffs_cleaned, fluxpath_name, diff) #convert dataframe to wide
fluxes_all_diffs_cleaned_wide <- column_to_rownames(fluxes_all_diffs_cleaned_wide, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names
#qvals
fluxes_qval_diffs_cleaned <- subset(fluxes_qval_diffs, select = -c(...1, Difference_direction, normal, tumor, fluxpath_description) )
#oh FFS, have to add underscores to flux descriptions or they won't get graphed. 
fluxes_qval_diffs_cleaned$fluxpath_name <- gsub(" ", "_", fluxes_qval_diffs_cleaned$fluxpath_name)
fluxes_qval_diffs_cleaned_wide <- spread(fluxes_qval_diffs_cleaned, fluxpath_name, diff) #convert dataframe to wide
fluxes_qval_diffs_cleaned_wide <- column_to_rownames(fluxes_qval_diffs_cleaned_wide, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names
#pvals
fluxes_pval_diffs_cleaned <- subset(fluxes_pval_diffs, select = -c(...1, Difference_direction, normal, tumor, fluxpath_description) )
#oh FFS, have to add underscores to flux descriptions or they won't get graphed. 
fluxes_pval_diffs_cleaned$fluxpath_name <- gsub(" ", "_", fluxes_pval_diffs_cleaned$fluxpath_name)
fluxes_pval_diffs_cleaned_wide <- spread(fluxes_pval_diffs_cleaned, fluxpath_name, diff) #convert dataframe to wide
fluxes_pval_diffs_cleaned_wide <- column_to_rownames(fluxes_pval_diffs_cleaned_wide, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names

#clean up
rm(fluxes_all, fluxes_all_diffs, fluxes_all_diffs_cleaned, fluxes_pval_diffs, fluxes_pval_diffs_cleaned, fluxes_qval_diffs, fluxes_qval_diffs_cleaned, top_fluxes)

##############Gene expression##############
############Step 1: import files (NEW host_gene_expression file from 05.14.version- fuso samples only)
host_gene_expression <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Common_growth_medium_runs/Jan_May2021_MICOM_GR_analyses/Correlating_GE_with_GRs_Burns_data/11.12.21_log2foldchange_gene_expression_top_125_metabolic_DEGs_matrix_Fuso_samples_only.csv") #the row_names note makes sure the Patient_blind_ids are set as row names
#note here that we've already filtered by the top 125 genes in the python script that generated this file; 
#see 05.14.21_calculating_log2foldchange_Burns_RNAseq_filtered_add_metabolism_filter_125_genes.py for that code

unique_PBIDs_host <- unique(host_gene_expression$Patient_Blind_ID)

#also, no need to do any dataframe manipulation. It's all ready to go!

##############CORRELATING##############\
#create a special DF for special qvalue; remove the PBID that's missing from q-value subset of flux data
host_gene_expression2 <- subset(host_gene_expression, Patient_Blind_ID != missing)
#set patient_blind_ids as rownames
host_gene_expression <- column_to_rownames(host_gene_expression, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names
host_gene_expression2 <- column_to_rownames(host_gene_expression2, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names

#host_gene_expression <- column_to_rownames(host_gene_expression, var= "Patient_Blind_ID") #already done

#order dataframes by row name (Patient_Blind_ID)
host_gene_expression <- host_gene_expression[ order(as.numeric(row.names(host_gene_expression))), ]
host_gene_expression2 <- host_gene_expression2[ order(as.numeric(row.names(host_gene_expression2))), ]

#fluxes_top_diffs_cleaned <- fluxes_top_diffs_cleaned[ order(as.numeric(row.names(fluxes_top_diffs_cleaned))), ]

#crud. Host gene expression and fluxes have different sizes. why....
#Fixed- had to redo all the filtering 

#Doing Spearman correlations with Benjamini-Hochberg corrections
#setting CI to true; makes things take (a second or two) longer to run, but CI output is nice.
fluxpath_correlations_all <- corr.test(fluxes_all_diffs_cleaned_wide, host_gene_expression, use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5)
fluxpath_correlations_qvals <- corr.test(fluxes_qval_diffs_cleaned_wide, host_gene_expression2, use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5)
fluxpath_correlations_pvals <- corr.test(fluxes_pval_diffs_cleaned_wide, host_gene_expression, use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5)

#NOTE: if you don't get rid of patient_blind_ID as a column, it'll throw up an error: Error in cor(x, y, use = use, method = method) : 'x' must be numeric
#b/c patient_blind_id is set as a factor on import

################Filtering and graphing##############
#step 1: get correlations_round2$p into a longform dataframe with the columns Taxon, Metabolite, and P-value
fluxpath_correlations_all_corrected <- as.data.frame(fluxpath_correlations_all$p.adj) #this makes a dataframe
fluxpath_correlations_all_corrected <- tibble::rownames_to_column(fluxpath_correlations_all_corrected, "fluxpath_name") #THANKS DPLYR, hero of R!
fluxpath_correlations_all_uncorrected <- as.data.frame(fluxpath_correlations_all$p) #this makes a dataframe
fluxpath_correlations_all_uncorrected <- tibble::rownames_to_column(fluxpath_correlations_all_uncorrected, "fluxpath_name") #THANKS DPLYR, hero of R!
#
fluxpath_correlations_pvals_corrected <- as.data.frame(fluxpath_correlations_pvals$p.adj) #this makes a dataframe
fluxpath_correlations_pvals_corrected <- tibble::rownames_to_column(fluxpath_correlations_pvals_corrected, "fluxpath_name") #THANKS DPLYR, hero of R!
fluxpath_correlations_pvals_uncorrected <- as.data.frame(fluxpath_correlations_pvals$p) #this makes a dataframe
fluxpath_correlations_pvals_uncorrected <- tibble::rownames_to_column(fluxpath_correlations_pvals_uncorrected, "fluxpath_name") #THANKS DPLYR, hero of R!
#
fluxpath_correlations_qvals_corrected <- as.data.frame(fluxpath_correlations_qvals$p.adj) #this makes a dataframe
fluxpath_correlations_qvals_corrected <- tibble::rownames_to_column(fluxpath_correlations_qvals_corrected, "fluxpath_name") #THANKS DPLYR, hero of R!
fluxpath_correlations_qvals_uncorrected <- as.data.frame(fluxpath_correlations_qvals$p) #this makes a dataframe
fluxpath_correlations_qvals_uncorrected <- tibble::rownames_to_column(fluxpath_correlations_qvals_uncorrected, "fluxpath_name") #THANKS DPLYR, hero of R!

################Filtering and graphing##############
#step 1: get correlations into a longform dataframe with the columns Taxon, Metabolite, and P-value
fluxpath_correlations_all_corrected_long <- fluxpath_correlations_all_corrected %>%
  melt(id.var="fluxpath_name") %>%
  arrange(fluxpath_name, variable) 
fluxpath_correlations_all_corrected_long <- rename(fluxpath_correlations_all_corrected_long, c("gene_id"="variable", "p_value_adj"="value")) #rename columns
#
fluxpath_correlations_all_uncorrected_long <- fluxpath_correlations_all_uncorrected %>%
  melt(id.var="fluxpath_name") %>%
  arrange(fluxpath_name, variable) 
fluxpath_correlations_all_uncorrected_long <- rename(fluxpath_correlations_all_uncorrected_long, c("gene_id"="variable", "p_value_raw"="value")) #rename columns
#
fluxpath_correlations_qvals_corrected_long <- fluxpath_correlations_qvals_corrected %>%
  melt(id.var="fluxpath_name") %>%
  arrange(fluxpath_name, variable) 
fluxpath_correlations_qvals_corrected_long <- rename(fluxpath_correlations_qvals_corrected_long, c("gene_id"="variable", "p_value_adj"="value")) #rename columns
#
fluxpath_correlations_qvals_uncorrected_long <- fluxpath_correlations_qvals_uncorrected %>%
  melt(id.var="fluxpath_name") %>%
  arrange(fluxpath_name, variable) 
fluxpath_correlations_qvals_uncorrected_long <- rename(fluxpath_correlations_qvals_uncorrected_long, c("gene_id"="variable", "p_value_raw"="value")) #rename columns
#
fluxpath_correlations_pvals_corrected_long <- fluxpath_correlations_pvals_corrected %>%
  melt(id.var="fluxpath_name") %>%
  arrange(fluxpath_name, variable) 
fluxpath_correlations_pvals_corrected_long <- rename(fluxpath_correlations_pvals_corrected_long, c("gene_id"="variable", "p_value_adj"="value")) #rename columns
#
fluxpath_correlations_pvals_uncorrected_long <- fluxpath_correlations_pvals_uncorrected %>%
  melt(id.var="fluxpath_name") %>%
  arrange(fluxpath_name, variable) 
fluxpath_correlations_pvals_uncorrected_long <- rename(fluxpath_correlations_pvals_uncorrected_long, c("gene_id"="variable", "p_value_raw"="value")) #rename columns

#step 2: merge dataframes
fluxpath_correlations_all_long_merged <- merge(fluxpath_correlations_all_corrected_long,
                                               fluxpath_correlations_all_uncorrected_long, by=c("fluxpath_name","gene_id"))
fluxpath_correlations_all_long_merged <- fluxpath_correlations_all_long_merged[order(fluxpath_correlations_all_long_merged$p_value_raw),] #sort by p-value
#
fluxpath_correlations_qvals_long_merged <- merge(fluxpath_correlations_qvals_corrected_long,
                                                     fluxpath_correlations_qvals_uncorrected_long, by=c("fluxpath_name","gene_id"))
fluxpath_correlations_qvals_long_merged <- fluxpath_correlations_qvals_long_merged[order(fluxpath_correlations_qvals_long_merged$p_value_raw),] #sort by p-value
#
fluxpath_correlations_pvals_long_merged <- merge(fluxpath_correlations_pvals_corrected_long,
                                                     fluxpath_correlations_pvals_uncorrected_long, by=c("fluxpath_name","gene_id"))
fluxpath_correlations_pvals_long_merged <- fluxpath_correlations_pvals_long_merged[order(fluxpath_correlations_pvals_long_merged$p_value_raw),] #sort by p-value
#save THESE to excel
list_of_datasets <- list("all_fluxpath_GE_corrs" = fluxpath_correlations_all_long_merged,
                         "qsig_fluxpath_GE_corrs" = fluxpath_correlations_qvals_long_merged,
                         "psig_fluxpath_GE_corrs" = fluxpath_correlations_pvals_long_merged)
write.xlsx(list_of_datasets, "12.21.21_Burns2015_metabolic_gene_expression_Fuso_only_fluxpath_correlations.xlsx", row.names= FALSE, col.names=TRUE)

#step 4: collect data to graph.
fluxpaths_and_ges_combined <-cbind(fluxes_all_diffs_cleaned_wide, host_gene_expression) #NOTE: we can only do the cbind b/c we have the same rows in the same order.

#step 5: graphing
#for graphing, data comes from fluxpaths_and_ges_combined, and p-values come from fluxpath_correlations_all_uncorrected_long
#first, need to remove (e) from all column headers
colnames(fluxpaths_and_ges_combined)<-gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(fluxpaths_and_ges_combined)))
#I am aware that this doesn't specifically remove (e), but it's the only way to get rid of the ellipses AND it gets rid of e, so...

#cytosine flux and GPT
ggplot(fluxpaths_and_ges_combined, aes(x=EX_csn, y=GPT)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("Glutamic-Pyruvate Transaminase vs. Cytosine exchange flux,"), subtitle= expression("Spearman correlation p.raw = 0.000189")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="log2(fold change flux), mmmol/(gDW * h),\n Cytosine exchange") +
  scale_y_continuous(name = "Glutamic-Pyruvate Transaminase,\nlog2(fold change gene expression)") 

#histidine and Aspartoacylase
ggplot(fluxpaths_and_ges_combined, aes(x=EX_his_L, y=ASPA)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("Aspartoacylase vs. L-histidine exchange flux,"), subtitle= expression("Spearman correlation p.raw = 0.000213")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="log2(fold change flux), mmmol/(gDW * h),\n L-histidine exchange") +
  scale_y_continuous(name = "Aspartoacylase,\nlog2(fold change gene expression)") 

#cytosine flux and GPT
ggplot(fluxpaths_and_ges_combined, aes(x=EX_csn, y=CA2)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("Carbonic Anhydrase 2 vs. Cytosine exchange flux,"), subtitle= expression("Spearman correlation p.raw = 0.000222")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="log2(fold change flux), mmmol/(gDW * h),\n Cytosine exchange") +
  scale_y_continuous(name = "Carbonic Anhydrase 2,\nlog2(fold change gene expression)") 

