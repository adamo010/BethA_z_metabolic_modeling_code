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

#for the generation of the data used in this script, see
#05.14.21_calculating_log2foldchange_Burns_RNAseq_filtered_add_metabolism_filter_125_genes.py (for GEs)
#09.29.21_Burns2015_fluxpath_analysis_at_median_subsystem_level.R (for fluxes)- NOTE HERE that we're using median values moving forward, instead
#of average values. 

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/") 

##############Fluxpaths##############
top_fluxes <- read_csv("09.29.21_filtered_median_microbial_flux_subsystems_all_results_Burns2015_data.csv")
fluxes_all <-read_csv("09.29.21_filtered_median_microbial_flux_subsystems_all_data_Burns2015_data.csv")

############step 2: filter fluxes_all 
####now, the filtering step. here we decide how many correlations we want to do.
#have two options: use all (69), or use only those with pvals > 0.1 (9) or use only those with pvals < 0.05 (6). Create vectors for each
all_flux_names <- pull(top_fluxes,flux_ids) #make a vector from a column:thanks dplyr!
pval0.1_flux_names <- all_flux_names[1:9] #select top 9 elements from the vector
pval0.05_flux_names <- all_flux_names[1:6] #select top 6 elements from the vector

#use that to filter on the new column
fluxes_all_diffs <- filter(fluxes_all, fluxpath_subsystem %in% all_flux_names) #thanks again dplyr!
fluxes_pval0.1_diffs <- filter(fluxes_all, fluxpath_subsystem %in% pval0.1_flux_names) #thanks again dplyr!
fluxes_pval0.05_diffs <- filter(fluxes_all, fluxpath_subsystem %in% pval0.05_flux_names) #thanks again dplyr!

#set Patient_Blind_ID as factor.
fluxes_all_diffs$Patient_Blind_ID <- as.factor(fluxes_all_diffs$Patient_Blind_ID)
fluxes_pval0.1_diffs$Patient_Blind_ID <- as.factor(fluxes_pval0.1_diffs$Patient_Blind_ID)
fluxes_pval0.05_diffs$Patient_Blind_ID <- as.factor(fluxes_pval0.05_diffs$Patient_Blind_ID)

############step 3: split dataframes into imports (+ve) and exports (-ve)
#NOT NEEDED HERE

#############step 4: NON ESSENTIAL: take a quick look at the graphs:
#Boxplots!
ggplot(data=fluxes_all_diffs, aes(x=fluxpath_subsystem, y= tumor, fill=fluxpath_subsystem)) +
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
fluxes_all_diffs[,c("log2_fluxes_tumor")] <- log2(fluxes_all_diffs[,c("tumor")])
fluxes_all_diffs[,c("log2_fluxes_normal")] <- log2(fluxes_all_diffs[,c("normal")])
#NaNs produced: is this bad? A Problem? Who knows?
#reader, it was a problem. That was not solved by adding a constant 1.
#oh, it's that you can't take the log of a negative value. So, the way to get around this is to....
#standby. need to check original code where I collapsed by subsystems.
#FIXED: couldn't figure out how to account for zeroes (reverse Rx fluxes), so I took medians instead of means within a subsystem.
fluxes_pval0.1_diffs[,c("log2_fluxes_tumor")] <- log2(fluxes_pval0.1_diffs[,c("tumor")])
fluxes_pval0.1_diffs[,c("log2_fluxes_normal")] <- log2(fluxes_pval0.1_diffs[,c("normal")])
#
fluxes_pval0.05_diffs[,c("log2_fluxes_tumor")] <- log2(fluxes_pval0.05_diffs[,c("tumor")])
fluxes_pval0.05_diffs[,c("log2_fluxes_normal")] <- log2(fluxes_pval0.05_diffs[,c("normal")])

#add a differences column
fluxes_all_diffs$log2_flux_diffs <- (fluxes_all_diffs$log2_fluxes_tumor - fluxes_all_diffs$log2_fluxes_normal)
fluxes_pval0.1_diffs$log2_flux_diffs <- (fluxes_pval0.1_diffs$log2_fluxes_tumor - fluxes_pval0.1_diffs$log2_fluxes_normal)
fluxes_pval0.05_diffs$log2_flux_diffs <- (fluxes_pval0.05_diffs$log2_fluxes_tumor - fluxes_pval0.05_diffs$log2_fluxes_normal)

#another way to do this
#fluxes_top_diffs[,c("log2_flux_diffs_V3")] <- log2(fluxes_top_diffs$fluxes_tumor / fluxes_top_diffs$fluxes_normal)

############Step 6: drop cols, convert dataframes to wide, and set row names as Patient_Blind_ID for correlating
fluxes_all_diffs_cleaned <- subset(fluxes_all_diffs, select = -c(...1, diff, Difference_direction, normal, tumor, log2_fluxes_normal, log2_fluxes_tumor) )
#oh FFS, have to add underscores to flux subsystems or they won't get graphed. 
fluxes_all_diffs_cleaned$fluxpath_subsystem <- gsub(" ", "_", fluxes_all_diffs_cleaned$fluxpath_subsystem)
fluxes_all_diffs_cleaned <- spread(fluxes_all_diffs_cleaned, fluxpath_subsystem, log2_flux_diffs) #convert dataframe to wide
fluxes_all_diffs_cleaned <- column_to_rownames(fluxes_all_diffs_cleaned, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names
#
fluxes_pval0.1_diffs_cleaned <- subset(fluxes_pval0.1_diffs, select = -c(...1, diff, Difference_direction, normal, tumor, log2_fluxes_normal, log2_fluxes_tumor) )
#oh FFS, have to add underscores to flux subsystems or they won't get graphed. 
fluxes_pval0.1_diffs_cleaned$fluxpath_subsystem <- gsub(" ", "_", fluxes_pval0.1_diffs_cleaned$fluxpath_subsystem)
fluxes_pval0.1_diffs_cleaned <- spread(fluxes_pval0.1_diffs_cleaned, fluxpath_subsystem, log2_flux_diffs) #convert dataframe to wide
fluxes_pval0.1_diffs_cleaned <- column_to_rownames(fluxes_pval0.1_diffs_cleaned, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names
#
fluxes_pval0.05_diffs_cleaned <- subset(fluxes_pval0.05_diffs, select = -c(...1, diff, Difference_direction, normal, tumor, log2_fluxes_normal, log2_fluxes_tumor) )
#oh FFS, have to add underscores to flux subsystems or they won't get graphed. 
fluxes_pval0.05_diffs_cleaned$fluxpath_subsystem <- gsub(" ", "_", fluxes_pval0.05_diffs_cleaned$fluxpath_subsystem)
fluxes_pval0.05_diffs_cleaned <- spread(fluxes_pval0.05_diffs_cleaned, fluxpath_subsystem, log2_flux_diffs) #convert dataframe to wide
fluxes_pval0.05_diffs_cleaned <- column_to_rownames(fluxes_pval0.05_diffs_cleaned, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names

##############Gene expression##############
############Step 1: import files (NEW host_gene_expression file from 05.04.version)
host_gene_expression <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Common_growth_medium_runs/Jan_May2021_MICOM_GR_analyses/Correlating_GE_with_GRs_Burns_data/05.14.21_log2foldchange_gene_expression_top_125_metabolic_DEGs_matrix.csv") #the row_names note makes sure the Patient_blind_ids are set as row names
#note here that we've already filtered by the top 125 genes in the python script that generated this file; 
#see 05.14.21_calculating_log2foldchange_Burns_RNAseq_filtered_add_metabolism_filter_125_genes.py for that code

#also, no need to do any dataframe manipulation. It's all ready to go!

##############CORRELATING##############
#set patient_blind_ids as rownames
host_gene_expression <- column_to_rownames(host_gene_expression, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names
#fluxes_top_diffs_cleaned <- column_to_rownames(fluxes_top_diffs_cleaned, var= "Patient_Blind_ID") #already done

#order dataframes by row name (Patient_Blind_ID)
host_gene_expression <- host_gene_expression[ order(as.numeric(row.names(host_gene_expression))), ]
#fluxes_top_diffs_cleaned <- fluxes_top_diffs_cleaned[ order(as.numeric(row.names(fluxes_top_diffs_cleaned))), ]

#Doing Spearman correlations with Benjamini-Hochberg corrections
#setting CI to true; makes things take (a second or two) longer to run, but CI output is nice.
fluxpath_correlations_all_pvals <- corr.test(fluxes_all_diffs_cleaned, host_gene_expression, use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5)
fluxpath_correlations_pval0.1_pvals <- corr.test(fluxes_pval0.1_diffs_cleaned, host_gene_expression, use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5)
fluxpath_correlations_pval0.05_pvals <- corr.test(fluxes_pval0.05_diffs_cleaned, host_gene_expression, use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5)

#NOTE: if you don't get rid of patient_blind_ID as a column, it'll throw up an error: Error in cor(x, y, use = use, method = method) : 'x' must be numeric
#b/c patient_blind_id is set as a factor on import

################Filtering and graphing##############
#step 1: get correlations_round2$p into a longform dataframe with the columns Taxon, Metabolite, and P-value
fluxpath_correlations_all_pvals_corrected <- as.data.frame(fluxpath_correlations_all_pvals$p.adj) #this makes a dataframe
fluxpath_correlations_all_pvals_corrected <- tibble::rownames_to_column(fluxpath_correlations_all_pvals_corrected, "fluxpath_subsystem") #THANKS DPLYR, hero of R!
fluxpath_correlations_all_pvals_uncorrected <- as.data.frame(fluxpath_correlations_all_pvals$p) #this makes a dataframe
fluxpath_correlations_all_pvals_uncorrected <- tibble::rownames_to_column(fluxpath_correlations_all_pvals_uncorrected, "fluxpath_subsystem") #THANKS DPLYR, hero of R!
#
fluxpath_correlations_pval0.1_pvals_corrected <- as.data.frame(fluxpath_correlations_pval0.1_pvals$p.adj) #this makes a dataframe
fluxpath_correlations_pval0.1_pvals_corrected <- tibble::rownames_to_column(fluxpath_correlations_pval0.1_pvals_corrected, "fluxpath_subsystem") #THANKS DPLYR, hero of R!
fluxpath_correlations_pval0.1_pvals_uncorrected <- as.data.frame(fluxpath_correlations_pval0.1_pvals$p) #this makes a dataframe
fluxpath_correlations_pval0.1_pvals_uncorrected <- tibble::rownames_to_column(fluxpath_correlations_pval0.1_pvals_uncorrected, "fluxpath_subsystem") #THANKS DPLYR, hero of R!
#
fluxpath_correlations_pval0.05_pvals_corrected <- as.data.frame(fluxpath_correlations_pval0.05_pvals$p.adj) #this makes a dataframe
fluxpath_correlations_pval0.05_pvals_corrected <- tibble::rownames_to_column(fluxpath_correlations_pval0.05_pvals_corrected, "fluxpath_subsystem") #THANKS DPLYR, hero of R!
fluxpath_correlations_pval0.05_pvals_uncorrected <- as.data.frame(fluxpath_correlations_pval0.05_pvals$p) #this makes a dataframe
fluxpath_correlations_pval0.05_pvals_uncorrected <- tibble::rownames_to_column(fluxpath_correlations_pval0.05_pvals_uncorrected, "fluxpath_subsystem") #THANKS DPLYR, hero of R!

################Filtering and graphing##############
#step 1: get correlations into a longform dataframe with the columns Taxon, Metabolite, and P-value
fluxpath_correlations_all_pvals_corrected_long <- fluxpath_correlations_all_pvals_corrected %>%
  melt(id.var="fluxpath_subsystem") %>%
  arrange(fluxpath_subsystem, variable) 
fluxpath_correlations_all_pvals_corrected_long <- rename(fluxpath_correlations_all_pvals_corrected_long, c("gene_id"="variable", "p_value_adj"="value")) #rename columns
#
fluxpath_correlations_all_pvals_uncorrected_long <- fluxpath_correlations_all_pvals_uncorrected %>%
  melt(id.var="fluxpath_subsystem") %>%
  arrange(fluxpath_subsystem, variable) 
fluxpath_correlations_all_pvals_uncorrected_long <- rename(fluxpath_correlations_all_pvals_uncorrected_long, c("gene_id"="variable", "p_value_raw"="value")) #rename columns
#
fluxpath_correlations_pval0.1_pvals_corrected_long <- fluxpath_correlations_pval0.1_pvals_corrected %>%
  melt(id.var="fluxpath_subsystem") %>%
  arrange(fluxpath_subsystem, variable) 
fluxpath_correlations_pval0.1_pvals_corrected_long <- rename(fluxpath_correlations_pval0.1_pvals_corrected_long, c("gene_id"="variable", "p_value_adj"="value")) #rename columns
#
fluxpath_correlations_pval0.1_pvals_uncorrected_long <- fluxpath_correlations_pval0.1_pvals_uncorrected %>%
  melt(id.var="fluxpath_subsystem") %>%
  arrange(fluxpath_subsystem, variable) 
fluxpath_correlations_pval0.1_pvals_uncorrected_long <- rename(fluxpath_correlations_pval0.1_pvals_uncorrected_long, c("gene_id"="variable", "p_value_raw"="value")) #rename columns
#
fluxpath_correlations_pval0.05_pvals_corrected_long <- fluxpath_correlations_pval0.05_pvals_corrected %>%
  melt(id.var="fluxpath_subsystem") %>%
  arrange(fluxpath_subsystem, variable) 
fluxpath_correlations_pval0.05_pvals_corrected_long <- rename(fluxpath_correlations_pval0.05_pvals_corrected_long, c("gene_id"="variable", "p_value_adj"="value")) #rename columns
#
fluxpath_correlations_pval0.05_pvals_uncorrected_long <- fluxpath_correlations_pval0.05_pvals_uncorrected %>%
  melt(id.var="fluxpath_subsystem") %>%
  arrange(fluxpath_subsystem, variable) 
fluxpath_correlations_pval0.05_pvals_uncorrected_long <- rename(fluxpath_correlations_pval0.05_pvals_uncorrected_long, c("gene_id"="variable", "p_value_raw"="value")) #rename columns
#step 2: merge dataframes
fluxpath_correlations_all_pvals_long_merged <- merge(fluxpath_correlations_all_pvals_corrected_long,
                                               fluxpath_correlations_all_pvals_uncorrected_long, by=c("fluxpath_subsystem","gene_id"))
fluxpath_correlations_all_pvals_long_merged <- fluxpath_correlations_all_pvals_long_merged[order(fluxpath_correlations_all_pvals_long_merged$p_value_raw),] #sort by p-value
#
fluxpath_correlations_pval0.1_pvals_long_merged <- merge(fluxpath_correlations_pval0.1_pvals_corrected_long,
                                                     fluxpath_correlations_pval0.1_pvals_uncorrected_long, by=c("fluxpath_subsystem","gene_id"))
fluxpath_correlations_pval0.1_pvals_long_merged <- fluxpath_correlations_pval0.1_pvals_long_merged[order(fluxpath_correlations_pval0.1_pvals_long_merged$p_value_raw),] #sort by p-value
#
fluxpath_correlations_pval0.05_pvals_long_merged <- merge(fluxpath_correlations_pval0.05_pvals_corrected_long,
                                                     fluxpath_correlations_pval0.05_pvals_uncorrected_long, by=c("fluxpath_subsystem","gene_id"))
fluxpath_correlations_pval0.05_pvals_long_merged <- fluxpath_correlations_pval0.05_pvals_long_merged[order(fluxpath_correlations_pval0.05_pvals_long_merged$p_value_raw),] #sort by p-value
#save THESE to excel
list_of_datasets <- list("all_flux_subsystems_GE_corrs" = fluxpath_correlations_all_pvals_long_merged,
                         "psig0.1_subsystems_GE_corrs" = fluxpath_correlations_pval0.1_pvals_long_merged,
                         "psig0.05_subsystems_GE_corrs" = fluxpath_correlations_pval0.05_pvals_long_merged)
write.xlsx(list_of_datasets, "09.29.21_Burns2015_metabolic_gene_expression_microbe_fluxpath_correlations.xlsx", row.names= FALSE, col.names=TRUE)

##################NOT DOING ANYTHING ELSE HERE##############
#step 4: collect data to graph.
fluxpaths_and_ges_combined <-cbind(fluxes_all_diffs_cleaned, host_gene_expression) #NOTE: we can only do the cbind b/c we have the same rows in the same order.

#step 5: graphing
#for graphing, data comes from fluxpaths_and_ges_combined, and p-values come from fluxpath_sig_corr_p_values_uncorrected_long
ggplot(fluxpaths_and_ges_combined, aes(x=Glutamate_metabolism, y=PDE5A)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("PDE5A vs. Glutamate_metabolism subsystem flux,"), subtitle= expression("Spearman correlation p.adj = 0.00269")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="log2(fold change flux), mmmol/(gDW * h),\n Glutamate_metabolism subsystem") +
  scale_y_continuous(name = "PDE5A,\nlog2(fold change gene expression)") 

ggplot(fluxpaths_and_ges_combined, aes(x=Starch_and_sucrose_metabolism, y=PPP1R3C)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("PPP1R3C vs. Starch and sucrose metabolism subsystem flux,"), subtitle= expression("Spearman correlation p = 0.07785")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="log2(fold change flux), mmmol/(gDW * h),\n Starch and sucrose metabolism subsystem") +
  scale_y_continuous(name = "PPP1R3C,\nlog2(fold change gene expression)") 

ggplot(fluxpaths_and_ges_combined, aes(x=Fatty_acid_oxidation, y=CDC14A)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("CDC14A vs. Fatty_acid_oxidation subsystem flux,"), subtitle= expression("Spearman correlation p = 0.000441")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="log2(fold change flux), mmmol/(gDW * h),\n Fatty_acid_oxidation subsystem") +
  scale_y_continuous(name = "CDC14A,\nlog2(fold change gene expression)") 

ggplot(fluxpaths_and_ges_combined, aes(x=Fatty_acid_oxidation, y=PRKG2)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("PRKG2 vs. Fatty_acid_oxidation subsystem flux,"), subtitle= expression("Spearman correlation p = 0.000725")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="log2(fold change flux), mmmol/(gDW * h),\n Fatty_acid_oxidation subsystem") +
  scale_y_continuous(name = "PRKG2,\nlog2(fold change gene expression)") 

ggplot(fluxpaths_and_ges_combined, aes(x=Purine_synthesis, y=PYGM)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("PYGM vs. Purine_synthesis subsystem flux,"), subtitle= expression("Spearman correlation p = 0.000882")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="log2(fold change flux), mmmol/(gDW * h),\n Purine_synthesis subsystem") +
  scale_y_continuous(name = "PYGM,\nlog2(fold change gene expression)") 

ggplot(fluxpaths_and_ges_combined, aes(x=	Vitamin_B6_metabolism, y=TGDS)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("TGDS vs. Vitamin_B6_metabolism subsystem flux,"), subtitle= expression("Spearman correlation p = 0.000904")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="log2(fold change flux), mmmol/(gDW * h),\n Vitamin_B6_metabolism subsystem") +
  scale_y_continuous(name = "TGDS,\nlog2(fold change gene expression)") 
