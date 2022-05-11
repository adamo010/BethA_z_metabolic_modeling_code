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

#Hello! Today I am doing 16s analysis on all samples from the metabolic modeling study

#this is based on 04.26.22_running_PCA_on_all_metabolic_modeling_datasets.R
#I'm editing it because I think there was a subsetting issue. 

#######Step 1: set wd
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal/") 

#######Step 2: import files and clean up
#NOTE: these tables have a) had contaminants removed, and b) been abundance-filtered to remove any taxa below 0.1% in all samples
#HEREIN LIES OUR PROBLEM! These are not actually the tax tables that we used for the model building.
#thus begins the deviation from 04.26 version.
#NOTE: these aren't actually the count tab. they're just a list of OTUs. We'll get to the count tab later. 

Burns_otu_list <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/06.23.20_summary_of_otus_with_models.csv")
Burns_otu_list <- Burns_otu_list %>% select(c(OTU_ID))
Hale_otu_list <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Hale_2018_model_picking/08.05.21_summary_of_Hale_OTUs_with_models.csv")
Hale_otu_list <- Hale_otu_list %>% select(c(OTU_ID))
Niccolai_otu_list<- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Niccolai_2020_model_picking/10.22.21_summary_of_Niccolai_OTUs_with_models.csv")
Niccolai_otu_list <- Niccolai_otu_list %>% select(c(OTU_ID))

#one important note about phyloseq is that it requires absolute counts, rather than abundance data. We have that data, fortunately! (Just) need to filter.
Burns_abs_otu_table <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal/Burns_absolute_OTU_table.csv")
Burns_abs_otu_table <- Burns_abs_otu_table %>% rename("#OTU_ID"="OTU_ID")
Hale_abs_otu_table <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Hale2018/sOTU_read_counts_by_sample_for_MICOM.csv")
Hale_abs_otu_table <- Hale_abs_otu_table %>% select(-c(...1))
Hale_abs_otu_table <- Hale_abs_otu_table %>% rename("OTU_ID_reads"="OTU_ID")
Hale_abs_otu_table <- column_to_rownames(Hale_abs_otu_table, var = "OTU_ID") #so odd, Hale abs_otu_table needs to be transposed.  
Hale_abs_otu_table <- as.data.frame(t(as.matrix(Hale_abs_otu_table)))
Hale_abs_otu_table <- rownames_to_column(Hale_abs_otu_table, var = "OTU_ID") 
Niccolai_abs_otu_table <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Niccolai2020/Niccolai2020_OTU_read_counts_by_sample.csv")

#filter abs_otu_tables by OTU_IDs present in count_tab df. Works pretty clean. Thanks dplyr!
Burns_count_tab <- filter(Burns_abs_otu_table, OTU_ID %in% Burns_otu_list$OTU_ID)
Hale_count_tab <- filter(Hale_abs_otu_table, OTU_ID %in% Hale_otu_list$OTU_ID)
Niccolai_count_tab <- filter(Niccolai_abs_otu_table, OTU_ID %in% Niccolai_otu_list$OTU_ID)
#a good sanity check is to make sure that the number of ROWS in count_tab and count_tab2 are the same.
#the number of columns will differ because count_tab contains extra taxonomic information at the end. Class_level_sort, etc. 

#so, the good news is that things look okay. The bad news is that taxonomy reporting is inconsistent.
Burns_tax_tab <- Burns_count_tab %>% select(c(OTU_ID, taxonomy))
remove <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__") #clean up taxonomy.
for (delim in remove){
  Burns_tax_tab <- Burns_tax_tab %>%
    mutate_at("taxonomy", str_replace, delim, "")
}
rm(remove,delim)
Burns_tax_tab <- Burns_tax_tab %>% separate(taxonomy, c('kingdom', 'phylum', 'class', 'order','family', 'genus', 'species'), sep=";")
Burns_tax_tab <- data.frame(lapply(Burns_tax_tab, as.character), stringsAsFactors=FALSE)
#find all instances of [ or ] in dataframe and remove them. Have to go by column, I guess. 
Burns_tax_tab$class <- gsub("\\[|]","",Burns_tax_tab$class) #use \\ because square brackets are a special character. | means "or"
Burns_tax_tab$order <- gsub("\\[|]","",Burns_tax_tab$order)
Burns_tax_tab$family <- gsub("\\[|]","",Burns_tax_tab$family)
Burns_tax_tab$genus <- gsub("\\[|]","",Burns_tax_tab$genus)
Burns_tax_tab$species <- gsub("\\[|]","",Burns_tax_tab$species)
#okay, need Hale taxonomy.
Hale_tax_tab <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Hale2018/Hale_2018_OTUs_with_matched_taxonomy.csv")
Hale_tax_tab <- filter(Hale_tax_tab, OTU_ID %in% Hale_count_tab$OTU_ID)
Hale_tax_tab <- Hale_tax_tab %>% 
  select(c("OTU_ID", 'Kingdom', 'Phylum', 'Class', 'Order','Family', 'Genus')) 
Hale_tax_tab <- Hale_tax_tab %>%
  rename_with(tolower) 
Hale_tax_tab <- Hale_tax_tab %>%
  rename_at(c("otu_id"),.funs=toupper)
#hopefully Niccolai will be easier
Niccolai_tax_tab <- Niccolai_count_tab %>% select(c('OTU_ID', 'Domain', 'Phylum', 'Class', 'Order','Family', 'Genus'))
Niccolai_tax_tab <- Niccolai_tax_tab %>%
  rename_with(tolower) 
Niccolai_tax_tab <- Niccolai_tax_tab %>%
  rename_at(c("otu_id"),.funs=toupper)

#import metadata; get all metadata onto SampleID, Patient_Blind_ID, Description 
Burns_sample_info_tab <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Burns_2015_sample_metadata.csv")
Hale_sample_info_tab <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Hale2018/Hale_2018_metadata_for_MICOM_samples.csv")
Hale_sample_info_tab <- Hale_sample_info_tab %>% rename(c("sampleIDs_from_read_counts"="SampleID", "host_subject_id"="Patient_Blind_ID",
                                                          "normal_adjacent_or_tumor_tissue_specimen"="Description"))
Niccolai_sample_info_tab <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Niccolai2020/Niccolai2020_all_metadata.csv")
Niccolai_sample_info_tab <- Niccolai_sample_info_tab %>% rename(c("Sample_ID"="SampleID"))

#clean up and move on
rm(Burns_abs_otu_table, Hale_abs_otu_table, Niccolai_abs_otu_table, Burns_otu_list, Hale_otu_list, Niccolai_otu_list)

###########Step 3: Build some phyloseq objects
#notes from the tutorial:
#otu_table - Works on any numeric matrix. You must also specify if the species are rows or columns
#sample_data - Works on any data.frame. The rownames must match the sample names in the otu_table if you plan to combine them as a phyloseq-object
#tax_table - Works on any character matrix. The rownames must match the OTU names (taxa_names) of the otu_table 
#if you plan to combine it with a phyloseq-object.

#sample_info (metadata) tab: set rownames as sampleids (the column names for count_tab2 and tax_tab)
Burns_sample_info_tab <- na.omit(Burns_sample_info_tab) #weird, extra blank rows in here. Should now have 88 rows
Burns_sample_info_tab <- column_to_rownames(Burns_sample_info_tab, var = "SampleID") #set sampleid  as rownames
Hale_sample_info_tab <- na.omit(Hale_sample_info_tab)  #weird, extra blank rows in here. Should now have 92 rows
Hale_sample_info_tab <- column_to_rownames(Hale_sample_info_tab, var = "SampleID")
#Niccolai_sample_info_tab <- na.omit(Niccolai_sample_info_tab)  #don't run this for Niccolai, as the dataset is a bit sparse
Niccolai_sample_info_tab <- column_to_rownames(Niccolai_sample_info_tab, var = "SampleID")

#tax_table tab: set rownames as otu_IDs
Burns_tax_tab <- column_to_rownames(Burns_tax_tab, var="OTU_ID")
Hale_tax_tab <- column_to_rownames(Hale_tax_tab, var="OTU_ID")
Niccolai_tax_tab <- column_to_rownames(Niccolai_tax_tab, var="OTU_ID")

#otu_table/count_tab: set rownamses as otu_id
Burns_count_tab <- column_to_rownames(Burns_count_tab, var="OTU_ID")
Burns_count_tab <- Burns_count_tab %>% select(-c("taxonomy")) #some reason there's an extra column here
Hale_count_tab <- column_to_rownames(Hale_count_tab, var="OTU_ID")
Niccolai_count_tab <- column_to_rownames(Niccolai_count_tab, var="OTU_ID")
Niccolai_count_tab <- Niccolai_count_tab %>% select(-c("Domain", "Phylum", "Class", "Order", "Family", "Genus")) #some reason there's an extra column here

#make phyloseq objects
ps_Burns <- phyloseq(otu_table(as.matrix(Burns_count_tab), taxa_are_rows=TRUE), 
                     sample_data(Burns_sample_info_tab), 
                     tax_table(as.matrix(Burns_tax_tab)))
ps_Hale <- phyloseq(otu_table(as.matrix(Hale_count_tab), taxa_are_rows=TRUE), 
                    sample_data(Hale_sample_info_tab), 
                    tax_table(as.matrix(Hale_tax_tab)))
ps_Niccolai <- phyloseq(otu_table(as.matrix(Niccolai_count_tab), taxa_are_rows=TRUE), 
                        sample_data(Niccolai_sample_info_tab), 
                        tax_table(as.matrix(Niccolai_tax_tab)))

#############Step 4: plotting alpha diversity
#all right, let's plot alpha diversity
BurnsAlpha <- plot_richness(ps_Burns, x="Description", measures=c("Shannon", "Simpson", "Chao1"), color="Description",
                            title="Burns 2015 dataset")+
  geom_boxplot(color="black", aes(fill= Description)) +
  geom_point(color="black") +
  theme_bw() +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
  theme(legend.position='none', axis.title.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("tumor", "normal")), label = "p.signif") +
  theme(strip.background = element_rect(color="black", size=3, linetype="blank"), strip.placement = "outside",
        strip.text.x = element_text(size = 10, color = "black", face = "bold"), plot.title = element_text(hjust=0.5))

HaleAlpha <- plot_richness(ps_Hale, x="Description", measures=c("Shannon", "Simpson", "Chao1"), 
                           color="Description", title="Hale 2018 dataset")+
  geom_boxplot(color="black", aes(fill= Description)) +
  geom_point(color="black") +
  theme_bw() +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
  theme(legend.position='none', axis.title.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("tumor", "normal")), label = "p.signif") +
  theme(strip.background = element_rect(color="black", size=3, linetype="blank"), strip.placement = "outside",
        strip.text.x = element_text(size = 10, color = "black", face = "bold"), plot.title = element_text(hjust=0.5))

NiccolaiAlpha <- plot_richness(ps_Niccolai, x="Description", measures=c("Shannon", "Simpson", "Chao1"), color="Description",
                               title="Niccolai 2020 dataset")+
  geom_boxplot(color="black", aes(fill= Description)) +
  geom_point(color="black") +
  theme_bw() +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
  theme(legend.position='none', axis.title.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("tumor", "normal")), label = "p.signif") +
  theme(strip.background = element_rect(color="black", size=3, linetype="blank"), strip.placement = "outside",
        strip.text.x = element_text(size = 10, color = "black", face = "bold"), plot.title = element_text(hjust=0.5))

#combine graphs
alpha_diversity <- ggarrange(BurnsAlpha,HaleAlpha,NiccolaiAlpha, 
                             ncol = 2, nrow = 2)

alpha_diversity #600x800

#clean up a bit
rm(BurnsAlpha, HaleAlpha, NiccolaiAlpha)

###########Step 5: Beta Diversity
# Transform data to proportions as appropriate for Bray-Curtis distances
# We're just going to use CLR. See https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
library(microbiome)
#clr-transform
ps_Burns_clr <- microbiome::transform(ps_Burns, "clr")
ps_Hale_clr <- microbiome::transform(ps_Hale, "clr")
ps_Niccolai_clr <- microbiome::transform(ps_Niccolai, "clr")
#ordinate by RDA
ps_Burns_clr_ord <- phyloseq::ordinate(ps_Burns_clr, "RDA")
ps_Hale_clr_ord <- phyloseq::ordinate(ps_Hale_clr, "RDA")
ps_Niccolai_clr_ord <- phyloseq::ordinate(ps_Niccolai_clr, "RDA")

#now we need some stats.
#let's run permanova
#Generate distance matrix
Burns_clr_dist_matrix <- phyloseq::distance(ps_Burns_clr, method = "euclidean") 
Hale_clr_dist_matrix <- phyloseq::distance(ps_Hale_clr, method = "euclidean") 
Niccolai_clr_dist_matrix <- phyloseq::distance(ps_Niccolai_clr, method = "euclidean") 
#ADONIS test
vegan::adonis(Burns_clr_dist_matrix ~ phyloseq::sample_data(ps_Burns_clr)$Description) #hey, p=0.002, R2= 0.0229
vegan::adonis(Hale_clr_dist_matrix ~ phyloseq::sample_data(ps_Hale_clr)$Description) # p=0.999, R2 = 0.00688
vegan::adonis(Niccolai_clr_dist_matrix ~ phyloseq::sample_data(ps_Niccolai_clr)$Description) #p=0.004, R2=0.01948

#The test from ADONIS can be confounded by differences in dispersion (or spread)…so we want to check this as well.
Burns_dispr <- vegan::betadisper(Burns_clr_dist_matrix, phyloseq::sample_data(ps_Burns_clr)$Description)
permutest(Burns_dispr) #p=0.185, so there is no difference in dispersion
Hale_dispr <- vegan::betadisper(Hale_clr_dist_matrix, phyloseq::sample_data(ps_Hale_clr)$Description)
permutest(Hale_dispr) #p=0.266
Niccolai_dispr <- vegan::betadisper(Niccolai_clr_dist_matrix, phyloseq::sample_data(ps_Niccolai_clr)$Description)
permutest(Niccolai_dispr) #p=0.302

phyloseq::plot_ordination(ps_Burns, ps_Burns_clr_ord, color="Description") + 
  geom_point(size = 2) +
  #coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = Description), linetype = 2)

#Now, graph
#NOTE: to look at other pcs, add axes=c(3,4) as an argument to phyloseq::plot_ordination; default is 1,2
Burns_pca <- phyloseq::plot_ordination(ps_Burns, ps_Burns_clr_ord, type="samples", color="Description") + 
  geom_point(size = 2) +
  stat_ellipse(aes(group = Description), linetype = 2)+
  theme_bw() +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
  labs(title="Burns 2015 dataset", subtitle= "PERMANOVA R^2 = 0.0223, p= 0.002") +
  theme(plot.title = element_text(size=14, hjust = 0.5), plot.subtitle=element_text(hjust=0.5)) +
  theme(plot.subtitle = ggtext::element_markdown())

Hale_pca <- phyloseq::plot_ordination(ps_Hale, ps_Hale_clr_ord, type="samples", color="Description") + 
  geom_point(size = 2) +
  stat_ellipse(aes(group = Description), linetype = 2) +
  theme_bw() +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
  labs(title="Hale 2018 dataset", color="Description", subtitle= "PERMANOVA R^2 = 0.00688, p= 0.999") +
  theme(plot.title = element_text(size=14, hjust = 0.5), plot.subtitle=element_text(hjust=0.5)) +
  theme(plot.subtitle = ggtext::element_markdown())

Niccolai_pca <- phyloseq::plot_ordination(ps_Niccolai, ps_Niccolai_clr_ord, type="samples", color="Description") + 
  geom_point(size = 2) +
  stat_ellipse(aes(group = Description), linetype = 2) +
  theme_bw() +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
  labs(title="Niccolai 2020 dataset", subtitle= "PERMANOVA R^2 = 0.01948, p= 0.006") +
  theme(plot.title = element_text(size=14, hjust = 0.5), plot.subtitle=element_text(hjust=0.5)) +
  theme(plot.subtitle = ggtext::element_markdown())

#combine graphs
pcas <- ggarrange(Burns_pca,Hale_pca,Niccolai_pca, 
                  labels=c("B.", "C.", "D."),
                  ncol = 3, nrow = 1,
                  common.legend = TRUE, legend = "bottom")
pcas #add , legend.position='bottom' to themes to get this to look nice. 1300 x 500

pcas2 <- ggarrange(Burns_pca,Hale_pca,Niccolai_pca, 
                   labels=c("B.", "C.", "D."),
                   ncol = 1, nrow = 3)
pcas2

#################Step 6: DEseq2
library("DESeq2")

#NOTE here, that deseq2 requires a vst transform, not a clr transform, so that's what we're using
Burns_dds = phyloseq_to_deseq2(ps_Burns, ~factor(Patient_Blind_ID) + factor(Description))
Burns_dds = estimateSizeFactors(Burns_dds)
Burns_dds = estimateDispersions(Burns_dds)
Burns_vst = getVarianceStabilizedData(Burns_dds) #do the vst transform
ps_Burns_vst= ps_Burns #duplicate our phyloseq object just in case
otu_table(ps_Burns_vst) <- otu_table(Burns_vst, taxa_are_rows = TRUE) 
#now has the results of DESeq2 variance-stabilization of counts instead of the original counts
Burns_dds = DESeq(Burns_dds, test="Wald", fitType="parametric")
Burns_res <- results(Burns_dds)
Burns_res <- Burns_res[order(Burns_res$padj),] #re-order by p-value
head(Burns_res)
Burns_sigtab = Burns_res[which(Burns_res$padj < 0.05), ]
#import taxonomic info
Burns_sigtab = cbind(as(Burns_sigtab, "data.frame"), as(tax_table(ps_Burns_vst)[rownames(Burns_sigtab), ], "matrix"))
head(Burns_sigtab)
#MUCH BETTER

Hale_dds = phyloseq_to_deseq2(ps_Hale, ~factor(Patient_Blind_ID) + factor(Description))
#Hale_dds = estimateSizeFactors(Hale_dds)
#Hale_dds = estimateDispersions(Hale_dds)
Hale_dds = DESeq(Hale_dds, test="Wald", fitType="parametric", sfType = "poscounts")
Hale_vst = getVarianceStabilizedData(Hale_dds) #do the vst transform
ps_Hale_vst= ps_Hale #duplicate our phyloseq object just in case
otu_table(ps_Hale_vst) <- otu_table(Hale_vst, taxa_are_rows = TRUE) 
#now has the results of DESeq2 variance-stabilization of counts instead of the original counts
Hale_res <- results(Hale_dds)
Hale_res <- Hale_res[order(Hale_res$padj),] #re-order by p-value
head(Hale_res)
Hale_sigtab = Hale_res[which(Hale_res$padj < 0.05), ]
#import taxonomic info
Hale_sigtab = cbind(as(Hale_sigtab, "data.frame"), as(tax_table(ps_Hale_vst)[rownames(Hale_sigtab), ], "matrix"))
head(Hale_sigtab)

Niccolai_dds = phyloseq_to_deseq2(ps_Niccolai, ~factor(Patient_Blind_ID) + factor(Description))
#Niccolai_dds = estimateSizeFactors(Niccolai_dds)
#Niccolai_dds = estimateDispersions(Niccolai_dds)
Niccolai_dds = DESeq(Niccolai_dds, test="Wald", fitType="parametric", sfType = "poscounts")
Niccolai_vst = getVarianceStabilizedData(Niccolai_dds) #do the vst transform
ps_Niccolai_vst= ps_Niccolai #duplicate our phyloseq object just in case
otu_table(ps_Niccolai_vst) <- otu_table(Niccolai_vst, taxa_are_rows = TRUE) 
#now has the results of DESeq2 variance-stabilization of counts instead of the original counts
Niccolai_res <- results(Niccolai_dds)
Niccolai_res <- Niccolai_res[order(Niccolai_res$padj),] #re-order by p-value
head(Niccolai_res)
Niccolai_sigtab = Niccolai_res[which(Niccolai_res$padj < 0.05), ]
#import taxonomic info
Niccolai_sigtab = cbind(as(Niccolai_sigtab, "data.frame"), as(tax_table(ps_Niccolai_vst)[rownames(Niccolai_sigtab), ], "matrix"))
head(Niccolai_sigtab)

###################Step 7: graph taxa
#looks better than 04.26 version, but there are A LOT of diff abundant taxa in Hale and Niccolai. Too many to graph, really. 
require("ggrepel")
library(viridis)
theme_set(theme_bw())
ggplot(Burns_sigtab, aes(x=genus, y=log2FoldChange, color=genus)) + geom_point(size=5) + 
  geom_hline(yintercept = 0, color = "black", size = 1) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust=0.5, size=10), axis.text.y= element_text(size=10)) +
  scale_y_continuous(name = "log2FoldChange in taxon abundance\n normal(base) vs tumor", limits=c(-15,10)) +
  scale_x_discrete(name = "Genus", guide = guide_axis(n.dodge=2)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  labs(title = expression("Burns 2015, differentially abundant taxa"), size=14) +
  theme(plot.title = element_text(hjust=0.5), legend.position='none') +
  geom_label_repel(aes(label = signif(padj, digits=3)), color="black", box.padding= 0.5, size=3, xlim = c(NA, Inf))

ggplot(Hale_sigtab, aes(x=genus, y=log2FoldChange, color=genus)) + geom_point(size=5) + 
  geom_hline(yintercept = 0, color = "black", size = 1) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust=0.5, size=10), axis.text.y= element_text(size=10)) +
  scale_y_continuous(name = "log2FoldChange in taxon abundance\n normal(base) vs tumor", limits=c(-30,30)) +
  scale_x_discrete(name = "Genus", guide = guide_axis(n.dodge=2)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  labs(title = expression("Hale 2018, differentially abundant taxa"), size=14) +
  theme(plot.title = element_text(hjust=0.5), legend.position='none') +
  geom_label_repel(aes(label = signif(padj, digits=3)), color="black", box.padding= 0.5, size=3, xlim = c(NA, Inf)) 

ggplot(Niccolai_sigtab, aes(x=family, y=log2FoldChange, color=family)) + geom_point(size=5) + 
  geom_hline(yintercept = 0, color = "black", size = 1) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust=0.5, size=10), axis.text.y= element_text(size=10)) +
  scale_y_continuous(name = "log2FoldChange in taxon abundance\n normal(base) vs tumor", limits=c(-30,30)) +
  scale_x_discrete(name = "Family", guide = guide_axis(n.dodge=2)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  labs(title = expression("Niccolai 2020, differentially abundant taxa"), size=14) +
  theme(plot.title = element_text(hjust=0.5), legend.position='none') +
  geom_label_repel(aes(label = signif(padj, digits=3)), color="black", box.padding= 0.5, size=3, xlim = c(NA, Inf)) 

#can we just look at fuso? #subset(Burns_sigtab, phylum=="Firmicutes"
#Going to have to create a new dataframe where we can graph all these

Burns_fuso <- ggplot(subset(Burns_sigtab, phylum==" Fusobacteria"), aes(x=genus, y=log2FoldChange, color=genus)) + geom_point(size=5) + 
  geom_hline(yintercept = 0, color = "black", size = 1) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust=0.5, size=10), axis.text.y= element_text(size=10)) +
  scale_y_continuous(name = "log2FoldChange in taxon abundance\n normal(base) vs tumor", limits=c(-5,10)) +
  scale_x_discrete(name = "Genus", guide = guide_axis(n.dodge=2)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  labs(title = expression("Burns 2015, differentially abundant taxa"), size=14) +
  theme(plot.title = element_text(hjust=0.5), legend.position='none') +
  geom_label_repel(aes(label = signif(padj, digits=3)), color="black", box.padding= 0.5, size=3, xlim = c(NA, Inf)) +
  theme(axis.text.x = element_text(face = "italic"))

Hale_fuso <- ggplot(subset(Hale_sigtab, phylum == "Fusobacteria"), aes(x=genus, y=log2FoldChange, color=genus))+ geom_point(size=5) + 
  geom_hline(yintercept = 0, color = "black", size = 1) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust=0.5, size=10), axis.text.y= element_text(size=10)) +
  scale_y_continuous(name = "log2FoldChange in taxon abundance\n normal(base) vs tumor", limits=c(0,25)) +
  scale_x_discrete(name = "Genus") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  labs(title = expression("Hale 2018, Phylum Fusobacteria"), size=14) +
  theme(plot.title = element_text(hjust=0.5), legend.position='none') +
  geom_label_repel(aes(label = signif(padj, digits=3)), color="black", box.padding= 0.5, size=3, xlim = c(NA, Inf)) +
  theme(axis.text.x = element_text(face = "italic"))

Niccolai_fuso <- ggplot(subset(Niccolai_sigtab, phylum == "Fusobacteria"), aes(x=genus, y=log2FoldChange, color=genus))+ geom_point(size=5) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  theme(axis.text.x = element_text(hjust = 0.5, vjust=0.5, size=10), axis.text.y= element_text(size=10)) +
  scale_y_continuous(name = "log2FoldChange in taxon abundance\n normal(base) vs tumor", limits=c(0,10)) +
  scale_x_discrete(name = "Genus") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  labs(title = expression("Niccolai 2020, Phylum Fusobacteria"), size=14) +
  theme(plot.title = element_text(hjust=0.5), legend.position='none') +
  geom_label_repel(aes(label = signif(padj, digits=3)), color="black", box.padding= 0.5, size=3, xlim = c(NA, Inf)) +
  theme(axis.text.x = element_text(face = "italic"))

allfuso <- ggarrange(Burns_fuso,Hale_fuso,Niccolai_fuso, 
                   #labels=c(".", "C.", "D."),
                   ncol = 3, nrow = 1)
allfuso





