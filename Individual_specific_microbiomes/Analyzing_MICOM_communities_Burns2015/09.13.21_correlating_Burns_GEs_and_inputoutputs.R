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

#the goal here is to import differentially expressed genes and significantly different microbial input/outputs from the Burns RNAseq data 
#and my metabolic models (based on the Burns microbiome data) and correlate them.

#based on 05.14.21_correlating_Burns_metabolic_GEs_and_GRs_abundance_filtered_tumor_GRs_only_125_genes_20_taxa.R
#more complete analysese can be found there. 

#based on #and Date_TBA_correlating_Burns_125_metabolic_GEs_20_microbial_inputsoutputs.R

#for the generation of the data used in this script, see
#05.14.21_calculating_log2foldchange_Burns_RNAseq_filtered_add_metabolism_filter_125_genes.py (for GEs)
#09.10.21_Burns2015_inputoutput_analysis.R (for input/outputs)

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/") 

############step 1: import files
top_metabs <- read_csv("09.10.21_filtered_microbial_input_output_sig_results_Burns2015_data.csv")
metabs_all <-read_csv("09.10.21_filtered_microbial_input_output_data_Burns2015_data.csv")

#goal: filter this dataframe by the top_metabs dataframe. In that, metabolite_amount has the identifier metabolite_name

############step 2: filter metabs_all by top_metabs to pull out growth rates of top 10 taxa only
#take a vector of the metab_ids from top_metabs
top_metab_names <- pull(top_metabs,metab_ids) #make a vector from a column:thanks dplyr!
#use that to filter on the new column
metabs_top_diffs <- filter(metabs_all, metabolite_name %in% top_metab_names) #thanks again dplyr!

############step 3: split dataframes into imports (+ve) and exports (-ve)
#first, double check that tumor and normal columns always have the same sign
#create a new column that's the product of tumor and normal columns. If it's -ve anywhere, we're fucked.
metabs_top_diffs <- metabs_top_diffs %>% mutate(testo = normal * tumor)
#we're good. delete that sorting col.
metabs_top_diffs <- metabs_top_diffs[ , ! names(metabs_top_diffs) %in% c("testo")]
#well, if we're only using tumor values anyway, maybe it doesn't matter. We can either do differences, or just tumor
#first, split into two dataframes based on +ve and -ve values for tumor (i.e which metabolites are exported and imported)
imported_metabs_top_diffs <- metabs_top_diffs[which(metabs_top_diffs$tumor > 0),]
exported_metabs_top_diffs <- metabs_top_diffs[which(metabs_top_diffs$tumor < 0),]
#for exported metabolites, need to convert to +ve values.
exported_metabs_top_diffs <- exported_metabs_top_diffs %>% mutate(tumor= abs(tumor)) #absolute value
#this, conveniently, removes all NAN values also. Maybe we'll just run with this.

#clean up
rm(metabs_top_diffs)

#############step 4: NON ESSENTIAL: take a quick look at the graphs:
#Boxplots!
#ggplot(data=imported_metabs_top_diffs, aes(x=metabolite_name, y= tumor, fill=metabolite_name)) +
  #geom_boxplot() +
  #geom_jitter(position=position_jitter(0.2)) +
  #labs(title= "Imported metabolite fluxes, tumor samples") +
  #scale_x_discrete(name ="Metabolite name") +
  #theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  #scale_y_continuous(name = "Import flux, mmmol/(gDW * h)") +
  #theme(legend.position="none", plot.title = element_text(hjust=0.5, vjust=2, size=12))

#want to look at the graphs without the two big confounders, Adenosine monophosphate and mucin-type O-glycan No. 9
#imported_metabs_top_diffs_subset<-imported_metabs_top_diffs[!(imported_metabs_top_diffs$metabolite_name=="Adenosine monophosphate" | 
                                                                #imported_metabs_top_diffs$metabolite_name=="mucin-type O-glycan No 9" |
                                                                #imported_metabs_top_diffs$metabolite_name=="released mucin-type O-glycan No 9"),]

#ggplot(data=imported_metabs_top_diffs_subset, aes(x=metabolite_name, y= tumor, fill=metabolite_name)) +
  #geom_boxplot() +
  #geom_jitter(position=position_jitter(0.2)) +
  #labs(title= "Imported metabolite fluxes, tumor samples") +
  #scale_x_discrete(name ="Metabolite name") +
  #theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  #scale_y_continuous(name = "Import flux, mmmol/(gDW * h)") +
  #theme(legend.position="none", plot.title = element_text(hjust=0.5, vjust=2, size=12))

############step 5: calculate log2(import) and log2(export) fluxes for tumor samples only- needed b/c we're taking log2fold change of genes.
#calculate log2 values and put into new column
exported_metabs_top_diffs[,c("log2_metabs_tumor")] <- log2(exported_metabs_top_diffs[,c("tumor")])
imported_metabs_top_diffs[,c("log2_metabs_tumor")] <- log2(imported_metabs_top_diffs[,c("tumor")])

#clean up a couple of columns
exported_metabs_top_diffs <- subset(exported_metabs_top_diffs, select = -c(...1, sorting_col, diff, Difference_direction) )
imported_metabs_top_diffs <- subset(imported_metabs_top_diffs, select = -c(...1, sorting_col, diff, Difference_direction) )

############Step 6: drop cols, convert dataframes to wide, and set row names as Patient_Blind_ID for correlating
exported_metabs_top_diffs_cleaned <- subset(exported_metabs_top_diffs, select = -c(normal,tumor) )
#exported_metabs_top_diffs_cleaned <- spread(exported_metabs_top_diffs_cleaned, metabolite_name, log2_metabs_tumor) #convert dataframe to wide

imported_metabs_top_diffs_cleaned <- subset(imported_metabs_top_diffs, select = -c(normal,tumor) )
#imported_metabs_top_diffs_cleaned <- spread(imported_metabs_top_diffs_cleaned, metabolite_name, log2_metabs_tumor) #convert dataframe to wide

#cleanup
rm(exported_metabs_top_diffs, imported_metabs_top_diffs)

##############Gene expression##############

############Step 1: import files (NEW host_gene_expression file from 05.04.version)
host_gene_expression <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Common_growth_medium_runs/Jan_May2021_MICOM_GR_analyses/Correlating_GE_with_GRs_Burns_data/05.14.21_log2foldchange_gene_expression_top_125_metabolic_DEGs_matrix.csv") #the row_names note makes sure the Patient_blind_ids are set as row names
#note here that we've already filtered by the top 125 genes in the python script that generated this file; 
#see 05.14.21_calculating_log2foldchange_Burns_RNAseq_filtered_add_metabolism_filter_125_genes.py for that code

#also, no need to do any dataframe manipulation. It's all ready to go!

##############CORRELATING##############
#convert to wideform. Want metabolites as column names, rows as Patient_Blind_ids.
#exported_metabs_top_diffs_cleaned_wide <- spread(exported_metabs_top_diffs_cleaned, metabolite_name, log2_metabs_tumor)
#okay, normally the above would work, but there is only one exported metabolite in this dataset, so spread is confused.
#instead, JUST FOR THIS DF, drop metabolite_name column and rename log2_metabs_tumor to D-Galactose.
exported_metabs_top_diffs_cleaned_wide <- subset(exported_metabs_top_diffs_cleaned, select = -c(metabolite_name) )
exported_metabs_top_diffs_cleaned_wide <- rename(exported_metabs_top_diffs_cleaned_wide, "D-Galactose" = "log2_metabs_tumor")
#sort by Patient_Blind_ID, then set Patient_Blind_ID as rownames
exported_metabs_top_diffs_cleaned_wide <- exported_metabs_top_diffs_cleaned_wide[order(exported_metabs_top_diffs_cleaned_wide$Patient_Blind_ID),]
exported_metabs_top_diffs_cleaned_wide <- column_to_rownames(exported_metabs_top_diffs_cleaned_wide, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names

#this one works fine.
imported_metabs_top_diffs_cleaned_wide <- spread(imported_metabs_top_diffs_cleaned, metabolite_name, log2_metabs_tumor) 
#sort by Patient_Blind_ID, then set Patient_Blind_ID as rownames
imported_metabs_top_diffs_cleaned_wide <- imported_metabs_top_diffs_cleaned_wide[order(imported_metabs_top_diffs_cleaned_wide$Patient_Blind_ID),]
imported_metabs_top_diffs_cleaned_wide <- column_to_rownames(imported_metabs_top_diffs_cleaned_wide, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names

#exported_metabs_top_diffs_cleaned_wide2 <- exported_metabs_top_diffs_cleaned_wide[ order(row.names(exported_metabs_top_diffs_cleaned_wide)), ]
#imported_metabs_top_diffs_cleaned_wide <- imported_metabs_top_diffs_cleaned_wide[ order(as.numeric(row.names(imported_metabs_top_diffs_cleaned_wide))), ]

rm(exported_metabs_top_diffs_cleaned, exported_metabs_top_diffs_cleaned)

#SPECIAL NOTE for exported values- there are only five patient_blind_ids. So, need to filter host_gene_expression by those patient_blind_ids
#first, get the ids
export_row_names <- row.names(exported_metabs_top_diffs_cleaned_wide) 
#now, filter host_gene_expression by export_row_names
host_gene_expression_trunc <- subset(host_gene_expression, Patient_Blind_ID %in% export_row_names)
#fix up all host gene expression dfs for correlating
host_gene_expression <- column_to_rownames(host_gene_expression, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names
host_gene_expression <- host_gene_expression[ order(as.numeric(row.names(host_gene_expression))), ] #order dataframes by row name (Patient_Blind_ID)
host_gene_expression_trunc <- column_to_rownames(host_gene_expression_trunc, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names
host_gene_expression_trunc <- host_gene_expression_trunc[ order(as.numeric(row.names(host_gene_expression_trunc))), ] #order dataframes by row name (Patient_Blind_ID)

#extra step: HOPEFULLY a temporary one, since eventually we will have all the data. But need to remove a few Patient_Blind_IDs from host_gene_expression,
#as these are also currently missing from metabs_top_diffs_cleaned
#first, have to reset all patient_blind_ids from row names to columns (gross)
#host_gene_expression <- tibble::rownames_to_column(host_gene_expression, "Patient_Blind_ID") #bless dplyr
#exported_metabs_top_diffs_cleaned <- tibble::rownames_to_column(exported_metabs_top_diffs_cleaned, "Patient_Blind_ID")
#truncated_host_gene_expression <- host_gene_expression[host_gene_expression$Patient_Blind_ID 
                                                       #%in% exported_metabs_top_diffs_cleaned$Patient_Blind_ID,]
#exported_metabs_top_diffs_cleaned <- column_to_rownames(exported_metabs_top_diffs_cleaned, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names
#truncated_host_gene_expression2 <- truncated_host_gene_expression[,-1]
#rownames(truncated_host_gene_expression2) <- truncated_host_gene_expression[,1]
#a couple notes:
#don't need to do imported metabolites because they should have the same patient_blind_ids as exported metabolites
#also, had to do extra crap for gene expression. no idea why, kept getting a "Error: `.data` must be a data frame without row names." error otherwise.

#Doing Spearman correlations with Benjamini-Hochberg corrections
#setting CI to true; makes things take (a second or two) longer to run, but CI output is nice.
#YOU KNOW WHAT?!?? I'm not doing corrected p-values anymore. They're all 1. what do you want from me. See 05.14.21_correlating if you want it
#NOTE: if you don't get rid of patient_blind_ID as a column, it'll throw up an error: Error in cor(x, y, use = use, method = method) : 'x' must be numeric
#b/c patient_blind_id is set as a factor on import

import_correlations_uncorr_pvals <- corr.test(imported_metabs_top_diffs_cleaned_wide, host_gene_expression, use = "pairwise",method="spearman",adjust="none", alpha=.05,ci=TRUE,minlength=5)
export_correlations_uncorr_pvals <- corr.test(exported_metabs_top_diffs_cleaned_wide, host_gene_expression_trunc, use = "pairwise",method="spearman",adjust="none", alpha=.05,ci=TRUE,minlength=5)

#eggcellent.

#COMMENTED OUT: update date and save.
list_of_export_datasets <-
list("R_values" = export_correlations_uncorr_pvals$r, "P_values" = export_correlations_uncorr_pvals$p, "CIs" = export_correlations_uncorr_pvals$ci)
write.xlsx(list_of_export_datasets, "09.14.21_Burns2015_metabolic_gene_expression_microbe_exported_metabs_uncorrected_filtered.xlsx", row.names= TRUE, col.names=TRUE)
list_of_import_datasets <-
list("R_values" = import_correlations_uncorr_pvals$r, "P_values" = import_correlations_uncorr_pvals$p, "CIs" = import_correlations_uncorr_pvals$ci)
write.xlsx(list_of_import_datasets, "09.14.21_Burns2015_metabolic_gene_expression_microbe_imported_metabs_uncorrected_filtered.xlsx", row.names= TRUE, col.names=TRUE)

################Filtering and graphing##############
#step 1: get correlations_round2$p into a longform dataframe with the columns Taxon, Metabolite, and P-value
#corr_p_values_corrected <- as.data.frame(correlations_corrected_pvals$p) #this makes a dataframe
export_corr_p_values_uncorrected <- as.data.frame(export_correlations_uncorr_pvals$p) #this makes a dataframe
import_corr_p_values_uncorrected <- as.data.frame(import_correlations_uncorr_pvals$p) #this makes a dataframe

#Step 2: convert a matrix to a longform data frame
#first, set row_names as first column. I know, I know. After all that work to do the opposite...
export_corr_p_values_uncorrected <- tibble::rownames_to_column(export_corr_p_values_uncorrected, "taxon_id") #THANKS DPLYR, hero of R!
import_corr_p_values_uncorrected <- tibble::rownames_to_column(import_corr_p_values_uncorrected, "taxon_id") #THANKS DPLYR, hero of R!

#then, convert to longform
export_corr_p_values_uncorrected_long <- export_corr_p_values_uncorrected %>%
  melt(id.var="taxon_id") %>%
  arrange(taxon_id, variable) 
import_corr_p_values_uncorrected_long <- import_corr_p_values_uncorrected %>%
  melt(id.var="taxon_id") %>%
  arrange(taxon_id, variable) 

#Step 3: sort and filter by significant p-values
export_corr_p_values_uncorrected_long <- rename(export_corr_p_values_uncorrected_long, c("gene_id"="variable", "p_value"="value")) #rename columns
export_corr_p_values_uncorrected_long <- export_corr_p_values_uncorrected_long[order(export_corr_p_values_uncorrected_long$p_value),] #sort by p-value
export_sig_corr_p_values_uncorrected_long <- filter(export_corr_p_values_uncorrected_long, p_value < 0.050) #only include p-values less than 0.05. There are 86 total.
import_corr_p_values_uncorrected_long <- rename(import_corr_p_values_uncorrected_long, c("gene_id"="variable", "p_value"="value")) #rename columns
import_corr_p_values_uncorrected_long <- import_corr_p_values_uncorrected_long[order(import_corr_p_values_uncorrected_long$p_value),] #sort by p-value
import_sig_corr_p_values_uncorrected_long <- filter(import_corr_p_values_uncorrected_long, p_value < 0.050) #only include p-values less than 0.05. There are 86 total.
#save data
write.csv(as.data.frame(export_sig_corr_p_values_uncorrected_long), file="09.14.21_Burns2015_filtered_uncorr_exported_metabs_and_tumor_metabolic_GE_only_correlations_pvals.csv") #save this to csv
write.csv(as.data.frame(import_sig_corr_p_values_uncorrected_long), file="09.14.21_Burns2015_filtered_uncorr_imported_metabs_and_tumor_metabolic_GE_only_correlations_pvals.csv") #save this to csv

#step 4: collect data to graph.
exports_and_ges_combined <-cbind(exported_metabs_top_diffs_cleaned_wide, host_gene_expression_trunc) #NOTE: we can only do the cbind b/c we have the same rows in the same order.
imports_and_ges_combined <-cbind(imported_metabs_top_diffs_cleaned_wide, host_gene_expression) #NOTE: we can only do the cbind b/c we have the same rows in the same order.

#step 5a: graph exports
#for graphing, data comes from exports/imports_and_ges_combined, and p-values come from sig_corr_p_values_uncorrected_long
#chrissake, rename D-galactose as something that stupid R can track.
exports_and_ges_combined <- exports_and_ges_combined %>% rename(D_galactose = "D-Galactose") #horseshit

ggplot(exports_and_ges_combined, aes(x= D_galactose, y=AURKA)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("AURKA vs. D-Galactose,"), subtitle= expression("Spearman correlation p = 3.97^10-24")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="D-Galactose export,\nmmmol/(gDW * h)") +
  scale_y_continuous(name = "AURKA,\nlog2(fold change gene expression)") 

ggplot(exports_and_ges_combined, aes(x= D_galactose, y=PPAT)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("PPAT vs. D-Galactose,"), subtitle= expression("Spearman correlation p = 3.97^10-24")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="D-Galactose export,\nmmmol/(gDW * h)") +
  scale_y_continuous(name = "PPAT,\nlog2(fold change gene expression)") 

ggplot(exports_and_ges_combined, aes(x= D_galactose, y=PRMT3)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("PRMT3 vs. D-Galactose,"), subtitle= expression("Spearman correlation p = 3.97^10-24")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="D-Galactose export,\nmmmol/(gDW * h)") +
  scale_y_continuous(name = "PRMT3,\nlog2(fold change gene expression)") 

ggplot(exports_and_ges_combined, aes(x= D_galactose, y=BUB1B)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("BUB1B vs. D-Galactose,"), subtitle= expression("Spearman correlation p = 0.0374")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="D-Galactose export,\nmmmol/(gDW * h)") +
  scale_y_continuous(name = "BUB1B,\nlog2(fold change gene expression)") 

ggplot(exports_and_ges_combined, aes(x= D_galactose, y=CA2)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("CA2 vs. D-Galactose,"), subtitle= expression("Spearman correlation p = 0.0374")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="D-Galactose export,\nmmmol/(gDW * h)") +
  scale_y_continuous(name = "CA2,\nlog2(fold change gene expression)") 

ggplot(exports_and_ges_combined, aes(x= D_galactose, y=CDC25C)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("CDC25C vs. D-Galactose,"), subtitle= expression("Spearman correlation p = 0.0374")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="D-Galactose export,\nmmmol/(gDW * h)") +
  scale_y_continuous(name = "CDC25C,\nlog2(fold change gene expression)") 

#step 5b: graph imports
#double chrissake, rename ions as something that stupid R can track.
imports_and_ges_combined <- imports_and_ges_combined %>% rename(Fe2 = "Fe2+") #horseshit
imports_and_ges_combined <- imports_and_ges_combined %>% rename(adenosine_monophosphate = "Adenosine monophosphate") #horseshit


ggplot(imports_and_ges_combined, aes(x= Fe2, y=GALNT16)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("GALNT16 vs. Fe2+,"), subtitle= expression("Spearman correlation p = 0.000341")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="Fe2+ import,\nmmmol/(gDW * h)") +
  scale_y_continuous(name = "GALNT16,\nlog2(fold change gene expression)") 

ggplot(imports_and_ges_combined, aes(x= Fe2, y=GSTM5)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("GSTM5 vs. Fe2+,"), subtitle= expression("Spearman correlation p = 0.000475")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="Fe2+ import,\nmmmol/(gDW * h)") +
  scale_y_continuous(name = "GSTM5,\nlog2(fold change gene expression)") 

ggplot(imports_and_ges_combined, aes(x= Calcium, y=GSTP1)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("GSTP1 vs. Calcium,"), subtitle= expression("Spearman correlation p = 0.00288")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="Calcium import,\nmmmol/(gDW * h)") +
  scale_y_continuous(name = "GSTP1,\nlog2(fold change gene expression)") 

ggplot(imports_and_ges_combined, aes(x= Calcium, y=GSTM5)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("GSTM5 vs. Calcium,"), subtitle= expression("Spearman correlation p = 0.00920")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="Calcium import,\nmmmol/(gDW * h)") +
  scale_y_continuous(name = "GSTM5,\nlog2(fold change gene expression)") 

ggplot(imports_and_ges_combined, aes(x= Heme, y=PLPP3)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("PLPP3 vs. Heme,"), subtitle= expression("Spearman correlation p = 0.00983")) +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=12), plot.subtitle = element_text(hjust = 0.5, vjust =1, size=12)) +
  guides(fill="none") +
  scale_x_continuous(name ="Heme import,\nmmmol/(gDW * h)") +
  scale_y_continuous(name = "PLPP3,\nlog2(fold change gene expression)") 








