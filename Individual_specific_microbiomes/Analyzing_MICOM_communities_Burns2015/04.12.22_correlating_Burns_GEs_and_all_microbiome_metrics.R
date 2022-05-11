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

#based on 09.13.21_correlating_Burns_GEs_and_GRs.R and 09.13.21_correlating_Burns_GEs_and_fluxpaths.R

#for the generation of the data used in this script, see
#05.14.21_calculating_log2foldchange_Burns_RNAseq_filtered_add_metabolism_filter_125_genes.py (for GEs)
#09.10.21_Burns2015_indiv_spp_GR_analysis.R (for GRs)
#10.21.21_BUrns2015_fluxpath_analysis_exchanges_only.R (for fluxes)

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/") 

##############Growth rates##############
############step 1: import files
genus_stat_results <- read_csv("09.10.21_filtered_microbial_GR_sig_results_Burns2015_data.csv") #used to be top_genera
genus_stat_results <- select(genus_stat_results, -c(...1))
genus_growth_rates <-read_csv("09.10.21_filtered_microbial_GR_data_Burns2015_data.csv") #used to be growth_rates_all
genus_growth_rates <- select(genus_growth_rates, -c(...1))
#goal: filter this dataframe by the top_ten_genera dataframe. In that, genus_ids_point_one has the identifier family_genus
#remove all rows where Diff_direction is no change; this is because a) who cares, b) those GRs
#are zero anyway, and c) it gives a -inf value for log2GRs, which will screw with everything unless I remove it.
genus_growth_rates <- subset(genus_growth_rates, Difference_direction != "Zero change") 
#this removes ~40 rows, which is fine.

############step 2: filter growth_rates_all by top_ten_genera to pull out growth rates of top 20 taxa only
#take a vector of the genus_ids_point_one from top_ten_genera
top_genera_names <- pull(genus_stat_results,genus_ids_point_one) #make a vector from a column:thanks dplyr!
#use that to filter on the new column
growth_rates_top_genera <- filter(genus_growth_rates, famgen_col %in% top_genera_names) #thanks again dplyr!

############step 3: calculate log2(tumor v normal GRs)- needed b/c we're taking log2fold change of genes.
#add 1 to each column (needed for taking the log of things)
growth_rates_top_genera$growth_rate_normal_adj <- growth_rates_top_genera$mean_genus_GR_normal + 1
growth_rates_top_genera$growth_rate_tumor_adj <- growth_rates_top_genera$mean_genus_GR_tumor + 1
#calculate log2 values and put into new column
growth_rates_top_genera[,c("log2_growth_rate_normal","log2_growth_rate_tumor")] <- log2(growth_rates_top_genera[,c("growth_rate_normal_adj","growth_rate_tumor_adj")])

############Step 4: average across OTUs that appear more than once in the same Patient_Blind_ID
#note here that we're only doing tumor GRs, because that's all we care about at this point.
#first, create a new dataframe containing only relevant columns
growth_rates_subset <- growth_rates_top_genera[ , names(growth_rates_top_genera) %in% c("Patient_Blind_ID", "log2_growth_rate_tumor", "OTU_ID", "famgen_col")] 
#then create a new sorting column to get unique values for each famgen_col/Patient_Blind_ID combo
growth_rates_subset <- growth_rates_subset %>% 
  mutate(sorting_col = paste0(famgen_col, "_", Patient_Blind_ID)) 
#then average across famgen_col :
growth_rates_subset_averaged <- growth_rates_subset %>%        # Specify data frame
  group_by(sorting_col, Patient_Blind_ID, famgen_col) %>%    # Specify group indicator
  summarise_at(vars(log2_growth_rate_tumor),                   # Specify column
               list(av_log2GR_tumor = mean))         # Specify function
growth_rates_subset_averaged <- growth_rates_subset_averaged[ , ! names(growth_rates_subset_averaged) %in% c("sorting_col")] #drop sorting column

############Step 5: convert dataframe to wide and set row names as Patient_Blind_ID for correlating
growth_rates_subset_wide <- spread(growth_rates_subset_averaged, famgen_col, av_log2GR_tumor) #convert dataframe to wide
growth_rates_subset_wide <- column_to_rownames(growth_rates_subset_wide, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names
#hooray! Now we have a dataframe for correlating.

#clean up
rm(genus_growth_rates, growth_rates_subset, growth_rates_subset_averaged, growth_rates_top_genera, genus_stat_results, top_genera_names)

##############Fluxes##############
############step 1: import files
flux_values <- read_csv("10.25.21_filtered_sum_exchange_fluxes_only_alldata_Burns2015_data.csv") #used to be fluxes_top_diffs
flux_values <- select(flux_values, -c(...1))
flux_stats <-read_csv("10.25.21_filtered_sum_exchange_fluxes_only_stats_Burns2015_data.csv") #used to be fluxes_all
flux_stats <- select(flux_stats, -c(...1))

############step 2: filter by top 20 significantly different exchange pathways. This makes the most comparable results to GRs.
#note: here we're filtering by p-values, but the order is the same as if we filtered by p-value
flux_stats <- flux_stats %>% top_n(-20, by_flux_pvals)
#take a vector of the genus_ids_point_one from top_ten_genera
top_flux_names <- pull(flux_stats,flux_ids) #make a vector from a column:thanks dplyr!
#use that to filter on the new column
fluxes_top_pathways <- filter(flux_values, fluxpath_name %in% top_flux_names) #thanks again dplyr!
#then, rename columns
fluxes_top_pathways <- fluxes_top_pathways %>% rename(fluxes_normal = normal, fluxes_tumor = tumor) #rename columns
fluxes_top_pathways$Patient_Blind_ID <- as.factor(fluxes_top_pathways$Patient_Blind_ID)

############step 3: calculate log2(import) and log2(export) fluxes for tumor samples only- needed b/c we're taking log2fold change of genes.
#calculate log2 values and put into new column
#flux_values[,c("log2_fluxes_tumor")] <- log2(flux_values[,c("fluxes_tumor")])
#ugh, I forgot about the issue with negative fluxes. RECALL: exports are -ve, imports are +ve by default in MICOM
#need to add another column that marks each as such...
flux_values_nonegs <- fluxes_top_pathways  %>%
  mutate(tumor_flux_direction = case_when(
    fluxes_tumor >0 ~ "_import",
    fluxes_tumor < 0 ~ "_export"))
#then convert all negative fluxes to positive values.
flux_values_nonegs$fluxes_tumor <- abs(flux_values_nonegs$fluxes_tumor)
#NOW we can take the log value
flux_values_nonegs[,c("log2_fluxes_tumor")] <- log2(flux_values_nonegs[,c("fluxes_tumor")])
#finally, we want to make sure the import/export direction is preserved, so we'll append it to the fluxpath_name column
flux_values_nonegs$flux_ids_dirs <- str_c(flux_values_nonegs$fluxpath_name,"",flux_values_nonegs$tumor_flux_direction)
#haha if we use stringr vs tidyr, we don't have to rearrange the columns.

############Step 4: drop cols, convert dataframes to wide, and set row names as Patient_Blind_ID for correlating
flux_values_cleaned <- subset(flux_values_nonegs, select = -c(fluxpath_name,fluxes_normal,fluxes_tumor,diff,Difference_direction,tumor_flux_direction))
flux_values_wide <- spread(flux_values_cleaned, flux_ids_dirs, log2_fluxes_tumor) #convert dataframe to wide
flux_values_wide <- column_to_rownames(flux_values_wide, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names

#ugh..... we should only have 20 metabolite exchanges, but because of the import/export thing, we get 31. How to contend with the fact that different samples
#have net export vs net import....
#do we maybe just... make the log values negative? No, we can't do that. negative log2 values happen when the original value is less than one. 
#maybe I just... don't transform. 
flux_values_notfn_cleaned <- subset(fluxes_top_pathways, select = -c(fluxes_normal,diff,Difference_direction))
flux_values_notfn_wide <- spread(flux_values_notfn_cleaned, fluxpath_name, fluxes_tumor) #convert dataframe to wide
flux_values_notfn_wide <- column_to_rownames(flux_values_notfn_wide, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names

#clean up
rm(flux_stats, flux_values, top_flux_names, flux_values_cleaned, flux_values_nonegs, fluxes_top_pathways, flux_values_notfn_cleaned)

##############Gene expression##############
############Step 1: import files (NEW host_gene_expression file from 05.14.version- fuso samples only)
host_gene_expression <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Common_growth_medium_runs/Jan_May2021_MICOM_GR_analyses/Correlating_GE_with_GRs_Burns_data/05.14.21_log2foldchange_gene_expression_top_125_metabolic_DEGs_matrix.csv") #the row_names note makes sure the Patient_blind_ids are set as row names
host_gene_expression <- column_to_rownames(host_gene_expression, var = "Patient_Blind_ID") #set Patient_Blind_ID as row names
#note here that we've already filtered by the top 125 genes in the python script that generated this file; 
#see 05.14.21_calculating_log2foldchange_Burns_RNAseq_filtered_add_metabolism_filter_125_genes.py for that code

#also, no need to do any dataframe manipulation. It's all ready to go!

##############CORRELATING##############
##########Step 1: prep dataframes
#order dataframes by row name (Patient_Blind_ID)
host_gene_expression <- host_gene_expression[ order(as.numeric(row.names(host_gene_expression))), ]
growth_rates_subset_wide <- growth_rates_subset_wide[ order(as.numeric(row.names(growth_rates_subset_wide))), ]
#flux_values_wide not used
flux_values_notfn_wide <- flux_values_notfn_wide [ order(as.numeric(row.names(flux_values_notfn_wide ))), ]

##########Step 2: Spearman correlations with Benjamini-Hochberg corrections
#setting CI to true; makes things take (a second or two) longer to run, but CI output is nice.
#gr_correlations
gr_ge_correlations <- corr.test(growth_rates_subset_wide, host_gene_expression, use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5)
#fluxpath_correlations
fluxpath_ge_correlations <- corr.test(flux_values_notfn_wide, host_gene_expression, use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5)

###########Step 3: Filtering and graphing
#gr_correlations: corrected and uncorrected p-values, and r-values
gr_correlations_corrected <- as.data.frame(gr_ge_correlations$p.adj) #this makes a dataframe
gr_correlations_corrected <- tibble::rownames_to_column(gr_correlations_corrected, "taxon_name") #THANKS DPLYR, hero of R!
gr_correlations_uncorrected <- as.data.frame(gr_ge_correlations$p) #this makes a dataframe
gr_correlations_uncorrected <- tibble::rownames_to_column(gr_correlations_uncorrected, "taxon_name") #THANKS DPLYR, hero of R!
gr_correlations_rvals <- as.data.frame(gr_ge_correlations$r) #this makes a dataframe
gr_correlations_rvals <- tibble::rownames_to_column(gr_correlations_rvals, "taxon_name") #THANKS DPLYR, hero of R!

#fluxpath_correlations: corrected and uncorrected p-values, and r-values
fluxpath_correlations_corrected <- as.data.frame(fluxpath_ge_correlations$p.adj) #this makes a dataframe
fluxpath_correlations_corrected <- tibble::rownames_to_column(fluxpath_correlations_corrected, "fluxpath_name") #THANKS DPLYR, hero of R!
fluxpath_correlations_uncorrected <- as.data.frame(fluxpath_ge_correlations$p) #this makes a dataframe
fluxpath_correlations_uncorrected <- tibble::rownames_to_column(fluxpath_correlations_uncorrected, "fluxpath_name") #THANKS DPLYR, hero of R!
fluxpath_correlations_rvals <- as.data.frame(fluxpath_ge_correlations$r) #this makes a dataframe
fluxpath_correlations_rvals <- tibble::rownames_to_column(fluxpath_correlations_rvals, "fluxpath_name") #THANKS DPLYR, hero of R!

#Get correlations into a longform dataframe with the columns Taxon, Metabolite, and P-value and R-value
gr_correlations_corrected_long <- gr_correlations_corrected %>%
  melt(id.var="taxon_name") %>%
  arrange(taxon_name, variable) 
gr_correlations_corrected_long <- rename(gr_correlations_corrected_long, c("gene_id"="variable", "p_value_adj"="value")) #rename columns
#
gr_correlations_uncorrected_long <- gr_correlations_uncorrected %>%
  melt(id.var="taxon_name") %>%
  arrange(taxon_name, variable) 
gr_correlations_uncorrected_long <- rename(gr_correlations_uncorrected_long, c("gene_id"="variable", "p_value_raw"="value")) #rename columns
#
gr_correlations_rvals_long <- gr_correlations_rvals %>%
  melt(id.var="taxon_name") %>%
  arrange(taxon_name, variable) 
gr_correlations_rvals_long <- rename(gr_correlations_rvals_long, c("gene_id"="variable", "R_value"="value")) #rename columns
###
fluxpath_correlations_corrected_long <- fluxpath_correlations_corrected %>%
  melt(id.var="fluxpath_name") %>%
  arrange(fluxpath_name, variable) 
fluxpath_correlations_corrected_long <- rename(fluxpath_correlations_corrected_long, c("gene_id"="variable", "p_value_adj"="value")) #rename columns
#
fluxpath_correlations_uncorrected_long <- fluxpath_correlations_uncorrected %>%
  melt(id.var="fluxpath_name") %>%
  arrange(fluxpath_name, variable) 
fluxpath_correlations_uncorrected_long <- rename(fluxpath_correlations_uncorrected_long, c("gene_id"="variable", "p_value_raw"="value")) #rename columns
#
fluxpath_correlations_rvals_long <- fluxpath_correlations_rvals %>%
  melt(id.var="fluxpath_name") %>%
  arrange(fluxpath_name, variable) 
fluxpath_correlations_rvals_long <- rename(fluxpath_correlations_rvals_long, c("gene_id"="variable", "R_value"="value")) #rename columns

#####merge and save dataframes
gr_correlations_long_merged <- merge(merge(gr_correlations_corrected_long,
                                     gr_correlations_uncorrected_long, all=TRUE),
                                     gr_correlations_rvals_long, all=TRUE)
gr_correlations_long_merged <- gr_correlations_long_merged[order(gr_correlations_long_merged$p_value_raw),] #sort by p-value
fluxpath_correlations_long_merged <- merge(merge(fluxpath_correlations_corrected_long,
                                     fluxpath_correlations_uncorrected_long, all=TRUE),
                                     fluxpath_correlations_rvals_long, all=TRUE)
fluxpath_correlations_long_merged <- fluxpath_correlations_long_merged[order(fluxpath_correlations_long_merged$p_value_raw),] #sort by p-value

#########save the files, just in case
write.csv(as.data.frame(gr_correlations_long_merged), file="04.12.22_Burns2015_metabolic_gene_expression_by_gr_correlations.csv")
write.csv(as.data.frame(fluxpath_correlations_long_merged), file="04.12.22_Burns2015_metabolic_gene_expression_by_fluxpath_correlations.csv")

#clean up
rm(flux_values_wide, fluxpath_correlations_corrected, fluxpath_correlations_corrected_long, fluxpath_correlations_rvals, fluxpath_correlations_rvals_long)
rm(fluxpath_correlations_uncorrected, fluxpath_correlations_uncorrected_long, fluxpath_ge_correlations, gr_correlations_corrected, gr_correlations_corrected_long)
rm(gr_correlations_rvals, gr_correlations_rvals_long, gr_correlations_uncorrected, gr_correlations_uncorrected_long, gr_ge_correlations)

##############filter to only keep rows where unadjusted P<0.01
gr_correlations_long_merged_sigonly <- filter(gr_correlations_long_merged, p_value_raw < 0.01) #this gets us 20 observations; 0.05 gets us 107 correlations
fluxpath_correlations_long_merged_sigonly <- filter(fluxpath_correlations_long_merged, p_value_raw < 0.01) #this gets us to 14 observations; 0.05 gets us 106 correlations

#########save the files, just in case
write.csv(as.data.frame(gr_correlations_long_merged_sigonly), file="04.12.22_Burns2015_metabolic_gene_expression_by_gr_sigto0.01_correlations.csv")
write.csv(as.data.frame(fluxpath_correlations_long_merged_sigonly), file="04.12.22_Burns2015_metabolic_gene_expression_by_fluxpath_sigto0.01_correlations.csv")

#########collect data to graph.
grs_and_ges_combined <-cbind(growth_rates_subset_wide, host_gene_expression) #NOTE: we can only do the cbind b/c we have the same rows in the same order.
fluxpaths_and_ges_combined <-cbind(flux_values_notfn_wide, host_gene_expression) #NOTE: we can only do the cbind b/c we have the same rows in the same order.

######### graphing
#for graphing, data comes from fluxpaths_and_ges_combined, and p-values come from fluxpath_correlations_all_uncorrected_long
#first, need to remove (e) from all column headers
colnames(grs_and_ges_combined)<-gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(grs_and_ges_combined)))
colnames(fluxpaths_and_ges_combined)<-gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(fluxpaths_and_ges_combined)))
#I am aware that this doesn't specifically remove (e), but it's the only way to get rid of the ellipses AND it gets rid of e, so...

#proline flux and PRMT3
plot1 <- ggplot(fluxpaths_and_ges_combined, aes(x=EX_pro_L, y=PRMT3)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("Protein Arginine Methyltransferase 3\n expression vs. L-proline exchange flux,")) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, vjust=2, size=13), axis.text=element_text(size=10),
        axis.text.y = element_markdown(), plot.margin= margin(0.75, 0.25, 0.25, 0.25, "cm")) +
  guides(fill="none") +
  labs(x=expression(atop("Tumor L-proline exchange flux,", paste("mmmol/(gDW * h)"))),
       y=expression(atop("PRMT3", paste("log"[2]*"(fold change gene expression)")))) +
  annotate("text", x= 1000, y=0.4, label="Spearman correlation\n p.adj = 0.0362") +
  annotate("text", x= -1000, y=0.05, label="R=-0.604")

#proline flux and CYP19A1
plot2 <- ggplot(fluxpaths_and_ges_combined, aes(x=EX_pro_L, y=CYP19A1)) + 
  geom_point() +
  geom_smooth(method=lm) +
  labs(title = expression("Cytochrome P450 Family 19 Subfamily A Member 1\n expression vs. L-proline exchange flux")) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, vjust=1, size=13), axis.text=element_text(size=10),
        axis.text.y = element_markdown(), plot.margin = margin(0.75,0.25,0.25,0.25, "cm")) +
  guides(fill="none") +
  labs(x=expression(atop("Tumor L-proline exchange flux,", paste("mmmol/(gDW * h)"))),
       y=expression(atop("CYP19A1", paste("log"[2]*"(fold change gene expression)")))) +
  annotate("text", x=1000, y=0.4, label="Spearman correlation\n p.adj = 0.0583") +
  annotate("text", x= -1000, y=0.05, label="R=-0.574")

corr_fig <- ggarrange(plot1, plot2,
                      labels = c("C.", "D."),
                      ncol = 2, nrow = 1)
corr_fig

