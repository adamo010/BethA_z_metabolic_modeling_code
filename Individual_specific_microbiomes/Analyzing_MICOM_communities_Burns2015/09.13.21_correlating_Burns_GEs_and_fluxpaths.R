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

#for the generation of the data used in this script, see
#05.14.21_calculating_log2foldchange_Burns_RNAseq_filtered_add_metabolism_filter_125_genes.py (for GEs)
#09.10.21_Burns2015_fluxpath_analysis.R (for fluxes)

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/") 

##############Fluxpaths##############
############step 1: import files
top_fluxes <- read_csv("09.10.21_filtered_microbial_fluxes_top20_sig_results_Burns2015_data.csv")
fluxes_all <-read_csv("09.10.21_filtered_microbial_fluxpaths_top20_data_Burns2015_data.csv")
#goal: filter this dataframe by the top_metabs dataframe. In that, metabolite_amount has the identifier metabolite_name

############step 2: filter fluxes_all by top_fluxes to pull out fluxes of top 10 taxa only
#take a vector of the metab_ids from top_metabs
top_flux_names <- pull(top_fluxes,flux_ids) #make a vector from a column:thanks dplyr!
#use that to filter on the new column
fluxes_top_diffs <- filter(fluxes_all, fluxpath_name %in% top_flux_names) #thanks again dplyr!
#then, rename columns
fluxes_top_diffs <- fluxes_top_diffs %>% rename(fluxes_normal = normal, fluxes_tumor = tumor) #rename columns
fluxes_top_diffs$Patient_Blind_ID <- as.factor(fluxes_top_diffs$Patient_Blind_ID)

############step 3: split dataframes into imports (+ve) and exports (-ve)
#NOT NEEDED HERE

#############step 4: NON ESSENTIAL: take a quick look at the graphs:
#Boxplots!
ggplot(data=fluxes_top_diffs, aes(x=fluxpath_name, y= fluxes_tumor, fill=fluxpath_name)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  labs(title= "Flux path activity, tumor samples") +
  scale_x_discrete(name ="Fluxpath name") +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  scale_y_continuous(name = "Flux, mmmol/(gDW * h)") +
  theme(legend.position="none", plot.title = element_text(hjust=0.5, vjust=2, size=12))

#huh. not bad.

############step 5: calculate log2(import) and log2(export) fluxes for tumor samples only- needed b/c we're taking log2fold change of genes.
#calculate log2 values and put into new column
fluxes_top_diffs[,c("log2_fluxes_tumor")] <- log2(fluxes_top_diffs[,c("fluxes_tumor")])

############Step 6: drop cols, convert dataframes to wide, and set row names as Patient_Blind_ID for correlating
fluxes_top_diffs_cleaned <- subset(fluxes_top_diffs, select = -c(...1, fluxpath_description, sorting_col, diff, Difference_direction, fluxes_normal,fluxes_tumor) )
fluxes_top_diffs_cleaned <- spread(fluxes_top_diffs_cleaned, fluxpath_name, log2_fluxes_tumor) #convert dataframe to wide
fluxes_top_diffs_cleaned <- column_to_rownames(fluxes_top_diffs_cleaned, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names

##############Gene expression##############
############Step 1: import files (NEW host_gene_expression file from 05.04.version)
host_gene_expression <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Common_growth_medium_runs/Jan_May2021_MICOM_GR_analyses/Correlating_GE_with_GRs_Burns_data/05.14.21_log2foldchange_gene_expression_top_125_metabolic_DEGs_matrix.csv") #the row_names note makes sure the Patient_blind_ids are set as row names
#note here that we've already filtered by the top 125 genes in the python script that generated this file; 
#see 05.14.21_calculating_log2foldchange_Burns_RNAseq_filtered_add_metabolism_filter_125_genes.py for that code

#also, no need to do any dataframe manipulation. It's all ready to go!

##############CORRELATING##############
#set patient_blind_ids as rownames
host_gene_expression <- column_to_rownames(host_gene_expression, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names
fluxes_top_diffs_cleaned <- column_to_rownames(fluxes_top_diffs_cleaned, var= "Patient_Blind_ID")

#order dataframes by row name (Patient_Blind_ID)
host_gene_expression <- host_gene_expression[ order(as.numeric(row.names(host_gene_expression))), ]
fluxes_top_diffs_cleaned <- fluxes_top_diffs_cleaned[ order(as.numeric(row.names(fluxes_top_diffs_cleaned))), ]

#Doing Spearman correlations with Benjamini-Hochberg corrections
#setting CI to true; makes things take (a second or two) longer to run, but CI output is nice.
fluxpath_correlations_corr_pvals <- corr.test(fluxes_top_diffs_cleaned, host_gene_expression, use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5)
fluxpath_correlations_uncorr_pvals <- corr.test(fluxes_top_diffs_cleaned, host_gene_expression, use = "pairwise",method="spearman",adjust="none", alpha=.05,ci=TRUE,minlength=5)
#NOTE: if you don't get rid of patient_blind_ID as a column, it'll throw up an error: Error in cor(x, y, use = use, method = method) : 'x' must be numeric
#b/c patient_blind_id is set as a factor on import

################Filtering and graphing##############
#step 1: get correlations_round2$p into a longform dataframe with the columns Taxon, Metabolite, and P-value
#corr_p_values_corrected <- as.data.frame(correlations_corrected_pvals$p) #this makes a dataframe
fluxpath_corr_p_values_corrected <- as.data.frame(fluxpath_correlations_corr_pvals$p) #this makes a dataframe
fluxpath_corr_p_values_uncorrected <- as.data.frame(fluxpath_correlations_uncorr_pvals$p) #this makes a dataframe
#here, corr and uncorr makes no difference at all, so move forward with uncorr.

#Step 2: convert a matrix to a longform data frame
#first, set row_names as first column. I know, I know. After all that work to do the opposite...
fluxpath_corr_p_values_uncorrected <- tibble::rownames_to_column(fluxpath_corr_p_values_uncorrected, "fluxpath_id") #THANKS DPLYR, hero of R!

#then, convert to longform
fluxpath_corr_p_values_uncorrected_long <- fluxpath_corr_p_values_uncorrected %>%
  melt(id.var="fluxpath_id") %>%
  arrange(fluxpath_id, variable) 

#Step 3: sort and filter by significant p-values
fluxpath_corr_p_values_uncorrected_long <- rename(fluxpath_corr_p_values_uncorrected_long, c("gene_id"="variable", "p_value"="value")) #rename columns
fluxpath_corr_p_values_uncorrected_long <- fluxpath_corr_p_values_uncorrected_long[order(fluxpath_corr_p_values_uncorrected_long$p_value),] #sort by p-value
fluxpath_sig_corr_p_values_uncorrected_long <- filter(fluxpath_corr_p_values_uncorrected_long, p_value < 0.050) #only include p-values less than 0.05. There are 86 total.
write.csv(as.data.frame(fluxpath_sig_corr_p_values_uncorrected_long), file="09.14.21_Burns2015_filtered_uncorr_fluxpaths_and_tumor_metabolic_GE_only_correlations_pvals.csv") #save this to csv

#step 4: collect data to graph.
fluxpaths_and_ges_combined <-cbind(fluxes_top_diffs_cleaned, host_gene_expression) #NOTE: we can only do the cbind b/c we have the same rows in the same order.

#step 5: graphing
#for graphing, data comes from fluxpaths_and_ges_combined, and p-values come from fluxpath_sig_corr_p_values_uncorrected_long
ggplot(fluxpaths_and_ges_combined, aes(x=TRE6PS, y=UGP2)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("PLPP1 vs. TRE6PS,"), subtitle= expression("Spearman correlation p = 7.98E-5")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="TRE6PS flux,\nmmmol/(gDW * h)") +
  scale_y_continuous(name = "PLPP1,\nlog2(fold change gene expression)") 

ggplot(fluxpaths_and_ges_combined, aes(x=GMPS2, y=PPA1)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("PPA1 vs. GMPS2,"), subtitle= expression("Spearman correlation p = 0.00210")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="GMPS2 flux,\nmmmol/(gDW * h)") +
  scale_y_continuous(name = "PPA1,\nlog2(fold change gene expression)") 

ggplot(fluxpaths_and_ges_combined, aes(x=OXAte, y=PPA1)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("PPA1 vs. OXAte,"), subtitle= expression("Spearman correlation p = 0.00210")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="OXAte flux,\nmmmol/(gDW * h)") +
  scale_y_continuous(name = "PPA1,\nlog2(fold change gene expression)") 

ggplot(fluxpaths_and_ges_combined, aes(x=GMPS2, y=MTHFD2)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("MTHFD2 vs. GMPS2,"), subtitle= expression("Spearman correlation p = 0.00273")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="GMPS2 flux,\nmmmol/(gDW * h)") +
  scale_y_continuous(name = "MTHFD2,\nlog2(fold change gene expression)") 

ggplot(fluxpaths_and_ges_combined, aes(x=OXAte, y=MTHFD2)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("MTHFD2 vs. OXAte,"), subtitle= expression("Spearman correlation p = 0.00286")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="OXAte flux,\nmmmol/(gDW * h)") +
  scale_y_continuous(name = "MTHFD2,\nlog2(fold change gene expression)") 

ggplot(fluxpaths_and_ges_combined, aes(x=OXAte, y=MET)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("MET vs. OXAte,"), subtitle= expression("Spearman correlation p = 0.00344")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="OXAte flux,\nmmmol/(gDW * h)") +
  scale_y_continuous(name = "MET,\nlog2(fold change gene expression)") 




