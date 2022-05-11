library("dada2")
library("stringr")
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
library("phyloseq")
library("Biostrings")

#Today I am doing 16s analysis on Muehlbauer dataset from the metabolic modeling study

#this is based on 04.29.22_V2_running_PCA_on_all_metabolic_modeling_datasets.R

#######Step 1: set wd
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal/") 

#######Step 2: import files and clean up
Muehlbauer_otu_list <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Muehlbauer_2021_model_picking/12.20.21_summary_of_Muehlbauer_OTUs_with_models.csv")
Muehlbauer_otu_list <- Muehlbauer_otu_list %>% select(c(OTU_ID))

#one important note about phyloseq is that it requires absolute counts, rather than abundance data. We have that data, fortunately! (Just) need to filter.
Muehlbauer_abs_otu_table <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/ASV_tables_Muehlbauer2021/Muehlbauer2021_OTUs_with_matched_taxonomy.csv")
Muehlbauer_abs_otu_table <- Muehlbauer_abs_otu_table %>% rename(c("...1"="OTU_ID"))

#filter abs_otu_tables by OTU_IDs present in count_tab df. Works pretty clean. Thanks dplyr!
Muehlbauer_count_tab <- filter(Muehlbauer_abs_otu_table, OTU_ID %in% Muehlbauer_otu_list$OTU_ID)

#extract taxonomic information
Muehlbauer_tax_tab <- Muehlbauer_count_tab %>% select(c('OTU_ID', 'kingdom', 'phylum', 'class', 'order','family', 'genus', 'species'))

#clean up count tab
Muehlbauer_count_tab <- Muehlbauer_count_tab %>% select(-c('kingdom', 'phylum', 'class', 'order','family', 'genus', 'species', 'phylum_level_sort',
                                                           'class_level_sort','order_level_sort','family_level_sort', 'genus_level_sort',
                                                           'species_level_sort', 'taxonomy_orig'))
#import metadata
Muehlbauer_sample_info_tab <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/ASV_tables_Muehlbauer2021/Muehlbauer2021_metadata_16s_only.csv")
Muehlbauer_sample_info_tab <- Muehlbauer_sample_info_tab %>% rename(c("Sample_name"="SampleID", "Isolate"="Patient_Blind_ID",
                                                          "experimental_treatment"="Description"))

#clean up and move on
rm(Muehlbauer_abs_otu_table, Muehlbauer_otu_list)

###########Step 3: Build some phyloseq objects
#notes from the tutorial:
#otu_table - Works on any numeric matrix. You must also specify if the species are rows or columns
#sample_data - Works on any data.frame. The rownames must match the sample names in the otu_table if you plan to combine them as a phyloseq-object
#tax_table - Works on any character matrix. The rownames must match the OTU names (taxa_names) of the otu_table 
#if you plan to combine it with a phyloseq-object.

#sample_info (metadata) tab: set rownames as sampleids (the column names for count_tab2 and tax_tab)
Muehlbauer_sample_info_tab <- column_to_rownames(Muehlbauer_sample_info_tab, var = "SampleID") #set sampleid  as rownames

#tax_table tab: set rownames as otu_IDs
Muehlbauer_tax_tab <- column_to_rownames(Muehlbauer_tax_tab, var="OTU_ID")

#otu_table/count_tab: set rownamses as otu_id
Muehlbauer_count_tab <- column_to_rownames(Muehlbauer_count_tab, var="OTU_ID")

#make phyloseq objects
ps_Muehlbauer <- phyloseq(otu_table(as.matrix(Muehlbauer_count_tab), taxa_are_rows=TRUE), 
                     sample_data(Muehlbauer_sample_info_tab), 
                     tax_table(as.matrix(Muehlbauer_tax_tab)))

#############Step 4: plotting alpha diversity
#all right, let's plot alpha diversity
plot_richness(ps_Muehlbauer, x="Description", measures=c("Shannon", "Simpson", "Chao1"), color="Description",
                            title="Muehlbauer 2021 dataset")+
  geom_boxplot(color="black", aes(fill= Description)) +
  geom_point(color="black") +
  theme_bw() +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
  theme(legend.position='none', axis.title.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  #stat_compare_means(method = "wilcox.test", comparisons = list(c("tumor", "normal")), label = "p.signif") +
  theme(strip.background = element_rect(color="black", size=3, linetype="blank"), strip.placement = "outside",
        strip.text.x = element_text(size = 10, color = "black", face = "bold"), plot.title = element_text(hjust=0.5))

###########Step 5: Beta Diversity
# Transform data to proportions as appropriate for Bray-Curtis distances
# We're just going to use CLR. See https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
library(microbiome)
#clr-transform
ps_Muehlbauer_clr <- microbiome::transform(ps_Muehlbauer, "clr")
#ordinate by RDA
ps_Muehlbauer_clr_ord <- phyloseq::ordinate(ps_Muehlbauer_clr, "RDA")

#now we need some stats.
#let's run permanova
#Generate distance matrix
Muehlbauer_clr_dist_matrix <- phyloseq::distance(ps_Muehlbauer_clr, method = "euclidean") 
#ADONIS test
vegan::adonis(Muehlbauer_clr_dist_matrix ~ phyloseq::sample_data(ps_Muehlbauer_clr)$Description) #p=0.976
#The test from ADONIS can be confounded by differences in dispersion (or spread)…so we want to check this as well.
Muehlbauer_dispr <- vegan::betadisper(Muehlbauer_clr_dist_matrix, phyloseq::sample_data(ps_Muehlbauer_clr)$Description)
permutest(Muehlbauer_dispr) #p=0.732, so there is no difference in dispersion

#Now, graph
#NOTE: to look at other pcs, add axes=c(3,4) as an argument to phyloseq::plot_ordination; default is 1,2
Muehlbauer_pca <- phyloseq::plot_ordination(ps_Muehlbauer, ps_Muehlbauer_clr_ord, type="samples", color="Description") + 
  geom_point(size = 2) +
  stat_ellipse(aes(group = Description), linetype = 2)+
  theme_bw() +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
  labs(title="Muehlbauer 2021 dataset", subtitle= "PERMANOVA R^2 = 0.0119, p= 0.976") +
  theme(plot.title = element_text(size=14, hjust = 0.5), plot.subtitle=element_text(hjust=0.5)) +
  theme(plot.subtitle = ggtext::element_markdown())

#Well...... that's rough.
#let's rerun run permanova, by Patient_Blind_ID this time
#ADONIS test
vegan::adonis(Muehlbauer_clr_dist_matrix ~ phyloseq::sample_data(ps_Muehlbauer_clr)$Patient_Blind_ID) #aha, p=0.001, R2=0.972
#The test from ADONIS can be confounded by differences in dispersion (or spread)…so we want to check this as well.
Muehlbauer_dispr <- vegan::betadisper(Muehlbauer_clr_dist_matrix, phyloseq::sample_data(ps_Muehlbauer_clr)$Patient_Blind_ID)
permutest(Muehlbauer_dispr) #p=0.597, so there is no difference in dispersion

#Now, graph
#NOTE: to look at other pcs, add axes=c(3,4) as an argument to phyloseq::plot_ordination; default is 1,2
Muehlbauer_pca <- phyloseq::plot_ordination(ps_Muehlbauer, ps_Muehlbauer_clr_ord, type="samples", color="Patient_Blind_ID") + 
  geom_point(size = 2) +
  #stat_ellipse(aes(group = Patient_Blind_ID), linetype = 2)+
  theme_bw() +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
  labs(title="Muehlbauer 2021 dataset", subtitle= "PERMANOVA R^2 = 0.972 , p= 0.001") +
  theme(plot.title = element_text(size=14, hjust = 0.5), plot.subtitle=element_text(hjust=0.5)) +
  theme(plot.subtitle = ggtext::element_markdown())





