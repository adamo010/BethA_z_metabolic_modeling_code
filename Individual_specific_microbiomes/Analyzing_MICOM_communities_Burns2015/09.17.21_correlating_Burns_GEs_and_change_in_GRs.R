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

#the goal here is to import differentially expressed genes and significantly different microbial growth rates from the Burns RNAseq data and my metabolic models
#(based on the Burns microbiome data) and correlate them.

#the difference between this and 05.03.21_correlating_Burns_GEs_and_GRs.R is that here, we're not using log2foldchangeGRs; we're just using log2(tumor GR)
#the difference between this and 05.04.21 version is that we're using the top 250 most differentially expressed METABOLIC genes.
#the difference between this and 05.12 version is that we're using the top 125 most differentially expressed METABOLIC genes and the top 20 taxa

#loosely based on 04.06.21_correlating_RNAseq_and_GRs_Burns_data.R, but ONLY the correlating part
#and based on 04.23.21_correlating_RNAseq_and_GRs_Burns_data_V2.R
#and Date_TBA_correlating_Burns_125_metabolic_GEs_20_microbial_taxaGRs.R
#and 09.13.21_correlating_Burns_GEs_and_GRs.R

#for the generation of the data used in this script, see
#05.14.21_calculating_log2foldchange_Burns_RNAseq_filtered_add_metabolism_filter_125_genes.py (for GEs)
#09.10.21_Burns2015_indiv_spp_GR_analysis.R (for GRs)

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/") 

##############Growth rates##############
############step 1: import files
top_genera <- read_csv("09.10.21_filtered_microbial_GR_sig_results_Burns2015_data.csv")
growth_rates_all <-read_csv("09.10.21_filtered_microbial_GR_data_Burns2015_data.csv")
#goal: filter this dataframe by the top_ten_genera dataframe. In that, genus_ids_point_one has the identifier family_genus

#one thing I'm going to add here is removing all rows where Diff_direction is no change; this is because a) who cares, b) those GRs
#are zero anyway, and c) it gives a -inf value for log2GRs, which will screw with everything unless I remove it.
growth_rates_all <- subset(growth_rates_all, Difference_direction != "Zero change") 
#this removes ~40 rows, which is fine.

############step 2: filter growth_rates_all by top_ten_genera to pull out growth rates of top 10 taxa only
#take family and genus columns from growth_rates_all and merge them into a new column
#growth_rates_all <- growth_rates_all %>% 
  #mutate(family_genus = paste0(family, "_", genus)) #dplyr: use paste0 to create a new column, family_genus, from the values in the "family" column 
#plus "_" plus the values in the genus column of the dataframe growth_rates_all

##don't actually need above code;  kept famgen col from GR analysis R script

#take a vector of the genus_ids_point_one from top_ten_genera
top_genera_names <- pull(top_genera,genus_ids_point_one) #make a vector from a column:thanks dplyr!
#use that to filter on the new column
growth_rates_top_genera <- filter(growth_rates_all, famgen_col %in% top_genera_names) #thanks again dplyr!

############step 3: calculate log2(tumor v normal GRs)- needed b/c we're taking log2fold change of genes.
#add 1 to each column (needed for taking the log of things)
growth_rates_top_genera$growth_rate_normal_adj <- growth_rates_top_genera$mean_genus_GR_normal + 1
growth_rates_top_genera$growth_rate_tumor_adj <- growth_rates_top_genera$mean_genus_GR_tumor + 1
#calculate log2 values and put into new column
growth_rates_top_genera[,c("log2_growth_rate_normal","log2_growth_rate_tumor")] <- log2(growth_rates_top_genera[,c("growth_rate_normal_adj","growth_rate_tumor_adj")])

############Step 4: average across OTUs that appear more than once in the same Patient_Blind_ID
#NEW here: instead of log2TumorGRs, get the difference: this gives log2fold change.
growth_rates_top_genera$log2_gr_diffs <- (growth_rates_top_genera$log2_growth_rate_tumor - growth_rates_top_genera$log2_growth_rate_normal)

#first, create a new dataframe containing only relevant columns
growth_rates_subset <- growth_rates_top_genera[ , names(growth_rates_top_genera) %in% c("Patient_Blind_ID", "log2_gr_diffs", "OTU_ID", "famgen_col")] 
#then create a new sorting column to get unique values for each famgen_col/Patient_Blind_ID combo
growth_rates_subset <- growth_rates_subset %>% 
  mutate(sorting_col = paste0(famgen_col, "_", Patient_Blind_ID)) 
#then average across famgen_col :
growth_rates_subset_averaged <- growth_rates_subset %>%        # Specify data frame
  group_by(sorting_col, Patient_Blind_ID, famgen_col) %>%    # Specify group indicator
  summarise_at(vars(log2_gr_diffs),                   # Specify column
               list(av_log2_gr_diffs = mean))         # Specify function
growth_rates_subset_averaged <- growth_rates_subset_averaged[ , ! names(growth_rates_subset_averaged) %in% c("sorting_col")] #drop sorting column

############Step 5: convert dataframe to wide and set row names as Patient_Blind_ID for correlating
growth_rates_subset_wide <- spread(growth_rates_subset_averaged, famgen_col, av_log2_gr_diffs) #convert dataframe to wide
growth_rates_subset_wide <- column_to_rownames(growth_rates_subset_wide, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names
#hooray! Now we have a dataframe for correlating.

#clean up
rm(growth_rates_all, growth_rates_subset, growth_rates_subset_averaged, growth_rates_top_genera, top_genera, top_genera_names)

##############Gene expression##############

############Step 1: import files (NEW host_gene_expression file from 05.04.version)
host_gene_expression <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Common_growth_medium_runs/Jan_May2021_MICOM_GR_analyses/Correlating_GE_with_GRs_Burns_data/05.14.21_log2foldchange_gene_expression_top_125_metabolic_DEGs_matrix.csv") #the row_names note makes sure the Patient_blind_ids are set as row names
host_gene_expression <- column_to_rownames(host_gene_expression, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names
#note here that we've already filtered by the top 125 genes in the python script that generated this file; 
#see 05.14.21_calculating_log2foldchange_Burns_RNAseq_filtered_add_metabolism_filter_125_genes.py for that code

#also, no need to do any dataframe manipulation. It's all ready to go!

##############CORRELATING##############
#order dataframes by row name (Patient_Blind_ID)
host_gene_expression <- host_gene_expression[ order(as.numeric(row.names(host_gene_expression))), ]
growth_rates_subset_wide <- growth_rates_subset_wide[ order(as.numeric(row.names(growth_rates_subset_wide))), ]

#Doing Spearman correlations with Benjamini-Hochberg corrections
#setting CI to true; makes things take (a second or two) longer to run, but CI output is nice.
#correlations_corrected_pvals <- corr.test(growth_rates_subset_wide, host_gene_expression, use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5)

#Look at the output:
#correlations_corrected_pvals$r %>%
  #head() #prints R value output
#correlations_corrected_pvals$p %>%
  #head() #prints p-value output
#correlations_corrected_pvals$ci %>%
  #head() #prints nice output table of CIs and other things (except that all the names of things are cut off.)

#save the outputs
#list_of_datasets <-
  #list("R_values" = correlations_corrected_pvals$r, "P_values" = correlations_corrected_pvals$p, "CIs" = correlations_corrected_pvals$ci)
#write.xlsx(list_of_datasets, "09.13.21_Burns2015_metabolic_gene_expression_microbe_tumor_GR_correlations_corrected_filtered.xlsx", row.names= TRUE, col.names=TRUE) #edited pre-corr row order

#Rerun, with unadjusted p-values.
#weirdly, for this particular run, they're all the same. Whatever. Weird shit.
correlations_uncorr_pvals <- corr.test(growth_rates_subset_wide, host_gene_expression, use = "pairwise",method="spearman",adjust="none", alpha=.05,ci=TRUE,minlength=5)
correlations_uncorr_pvals$p %>%
  head() #prints p-value output

#list_of_datasets <-
  #list("R_values" = correlations_uncorr_pvals$r, "P_values" = correlations_uncorr_pvals$p, "CIs" = correlations_uncorr_pvals$ci)
#write.xlsx(list_of_datasets, "09.13.21_Burns2015_metabolic_gene_expression_microbe_tumor_GR_correlations_uncorrected_filtered.xlsx", row.names= TRUE, col.names=TRUE)

################Filtering and graphing##############
#step 1: get correlations_round2$p into a longform dataframe with the columns Taxon, Gene, and P-value
#corr_p_values_corrected <- as.data.frame(correlations_corrected_pvals$p) #this makes a dataframe
corr_p_values_uncorrected <- as.data.frame(correlations_uncorr_pvals$p) #this makes a dataframe

#Step 2: convert a matrix to a longform data frame
#first, set row_names as first column. I know, I know. After all that work to do the opposite...
#corr_p_values_corrected <- tibble::rownames_to_column(corr_p_values_corrected, "taxon_id") #THANKS DPLYR, hero of R!
corr_p_values_uncorrected <- tibble::rownames_to_column(corr_p_values_uncorrected, "taxon_id") #THANKS DPLYR, hero of R!
#then, convert to longform
#corr_p_values_corrected_long <- corr_p_values_corrected %>%
  #melt(id.var="taxon_id") %>%
  #arrange(taxon_id, variable) 

corr_p_values_uncorrected_long <- corr_p_values_uncorrected %>%
  melt(id.var="taxon_id") %>%
  arrange(taxon_id, variable) 

#Step 3: sort and filter by significant p-values
#corr_p_values_corrected_long <- rename(corr_p_values_corrected_long, c("gene_id"="variable", "p_value"="value")) #rename columns
#corr_p_values_corrected_long <- corr_p_values_corrected_long[order(corr_p_values_corrected_long$p_value),] #sort by p-value
#sig_corr_p_values_corrected_long <- filter(corr_p_values_corrected_long, p_value < 0.050) #only include p-values less than 0.05. There are 86 total.
#write.csv(as.data.frame(sig_corr_p_values_corrected_long ), file="09.13.21_Burns2015_filtered_multiplecorr_GR_and_tumor_metabolic_GE_only_correlations_pvalues.csv") #save this to csv

corr_p_values_uncorrected_long <- rename(corr_p_values_uncorrected_long, c("gene_id"="variable", "p_value"="value")) #rename columns
corr_p_values_uncorrected_long <- corr_p_values_uncorrected_long[order(corr_p_values_uncorrected_long$p_value),] #sort by p-value
sig_corr_p_values_uncorrected_long <- filter(corr_p_values_uncorrected_long, p_value < 0.050) #only include p-values less than 0.05. There are 86 total.
write.csv(as.data.frame(sig_corr_p_values_uncorrected_long ), file="09.17.21_Burns2015_filtered_uncorr_change_in_GR_and_tumor_metabolic_GE_only_correlations_pvalues.csv") #save this to csv

#step 4: collect data to graph.
grs_and_ges_combined <-cbind(growth_rates_subset_wide, host_gene_expression) #NOTE: we can only do the cbind b/c we have the same rows in the same order.

#step 5: graph.
#for graphing, data comes from grs_and_ges_combined, and p-values come from sig_corr_p_values_uncorrected_long

#GRAPH 1
ggplot(grs_and_ges_combined, aes(x=Oxalobacteraceae_Ralstonia, y=DPAGT1)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("DPAGT1 vs. "*italic("Ralstonia")), subtitle= expression("Spearman correlation p = 0.000907")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="Ralstonia,\nlog2(tumor-normal growth rate)") +
  scale_y_continuous(name = "DPAGT1,\nlog2(fold change gene expression)") 

#GRAPH 2
ggplot(grs_and_ges_combined, aes(x=Oxalobacteraceae_Ralstonia, y=ACP1)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("ACP1 vs. "*italic("Ralstonia")), subtitle= expression("Spearman correlation p = 0.00139")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="Ralstonia,\nlog2(tumor-normal growth rate)") +
  scale_y_continuous(name = "ACP1,\nlog2(fold change gene expression)") 

#GRAPH 3
ggplot(grs_and_ges_combined, aes(x=Desulfovibrionaceae_Desulfovibrio, y=ATR)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("ATR vs. "*italic("Desulfovibrio")), subtitle= expression("Spearman correlation p = 0.00224")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="Desulfovibrio,\nlog2(tumor-normal growth rate)") +
  scale_y_continuous(name = "ATR,\nlog2(fold change gene expression)") 

#GRAPH 4
ggplot(grs_and_ges_combined, aes(x=	Porphyromonadaceae_Porphyromonas, y=BUB1B)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("BUB1B vs. "*italic("Porphyromonas")), subtitle= expression("Spearman correlation p = 0.00329")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="Porphyromonas,\nlog2(tumor-normal growth rate)") +
  scale_y_continuous(name = "BUB1B,\nlog2(fold change gene expression)") 

#GRAPH 5
ggplot(grs_and_ges_combined, aes(x= Oxalobacteraceae_Ralstonia, y=FARSB)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("FARSB vs. "*italic("Ralstonia")), subtitle= expression("Spearman correlation p = 0.00342")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="Ralstonia,\nlog2(tumor-normal growth rate)") +
  scale_y_continuous(name = "FARSB,\nlog2(fold change gene expression)") 

#GRAPH 6
ggplot(grs_and_ges_combined, aes(x= Desulfovibrionaceae_Desulfovibrio, y=CDC25B)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("CDC25B vs. "*italic("Desulfovibrio")), subtitle= expression("Spearman correlation p = 0.00486")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="Desulfovibrio,\nlog2(tumor-normal growth rate)") +
  scale_y_continuous(name = "CDC25B,\nlog2(fold change gene expression)") 
