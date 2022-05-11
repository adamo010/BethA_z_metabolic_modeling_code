#Hello! Today I am doing 16s analysis on all samples from the metabolic modeling study

#this is based on 04.15.22_mississippi_CRC_microbiome_phyloseq_only.R

######step 1: load packages
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

rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal/") 

######step 2: import files and clean up
#NOTE: these tables have a) had contaminants removed, and b) been abundance-filtered to remove any taxa below 0.1% in all samples
Burns_count_tab <- read_csv("07.24.20_otu_table_relative_abundances.csv")
Burns_count_tab <- Burns_count_tab %>% select(-c(...1))
remove <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__") #clean up taxonomy.
for (delim in remove){
  Burns_count_tab <- Burns_count_tab %>%
    mutate_at("taxonomy_orig", str_replace, delim, "")
}
rm(remove,delim)
Burns_count_tab <- Burns_count_tab %>% rename("#OTU_ID"="OTU_ID")
Hale_count_tab <- read_csv("07.13.21_Hale2018_otu_table_abundance_filtered_0.1_any_sample_cutoff.csv")
Niccolai_count_tab <- read_csv("10.19.21_Niccolai2020_otu_table_abundance_filtered_0.1_any_sample_cutoff.csv")

#one important note about phyloseq is that it requires absolute counts, rather than abundance data. We have that data, fortunately! (Just) need to filter.
Burns_abs_otu_table <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal/Burns_absolute_OTU_table.csv")
Burns_abs_otu_table <- Burns_abs_otu_table %>% rename("#OTU_ID"="OTU_ID")
Hale_abs_otu_table <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Hale2018/sOTU_read_counts_by_sample_for_MICOM.csv")
Hale_abs_otu_table <- Hale_abs_otu_table %>% select(-c(...1))
Hale_abs_otu_table <- Hale_abs_otu_table %>% rename("OTU_ID_reads"="OTU_ID")
Hale_abs_otu_table <- column_to_rownames(Hale_abs_otu_table, var = "OTU_ID") #so odd, Hale abs_otu_table needs to be transposed.  
Hale_abs_otu_table <- as.data.frame(t(as.matrix(Hale_abs_otu_table)))
rownames(Hale_abs_otu_table) <- NULL #weirdly, have OTU_ID elsewhere in dataframe
Niccolai_abs_otu_table <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Niccolai2020/Niccolai2020_OTU_read_counts_by_sample.csv")

#filter abs_otu_tables by OTU_IDs present in count_tab df. Works pretty clean. Thanks dplyr!
Burns_count_tab2 <- filter(Burns_abs_otu_table, OTU_ID %in% Burns_count_tab$OTU_ID)
Hale_count_tab2 <- filter(Hale_abs_otu_table, OTU_ID %in% Hale_count_tab$OTU_ID)
Niccolai_count_tab2 <- filter(Niccolai_abs_otu_table, OTU_ID %in% Niccolai_count_tab$OTU_ID)
#a good sanity check is to make sure that the number of ROWS in count_tab and count_tab2 are the same.
#the number of columns will differ because count_tab contains extra taxonomic information at the end. Class_level_sort, etc. 

#now, we need to extract the taxonomy information from each dataframe. lucky for me, they're all organized in the same way. 
Burns_tax_tab <- Burns_count_tab %>% select(c(OTU_ID, taxonomy_orig))
Burns_tax_tab <- Burns_tax_tab %>% separate(taxonomy_orig, c('kingdom', 'phylum', 'class', 'order','family', 'genus', 'species'), sep=";")
Hale_tax_tab <- Hale_count_tab %>% select(c(OTU_ID, taxonomy_orig))
Hale_tax_tab <- Hale_tax_tab %>% separate(taxonomy_orig, c('kingdom', 'phylum', 'class', 'order','family', 'genus', 'species'), sep=";")
Niccolai_tax_tab <- Niccolai_count_tab %>% select(c(OTU_ID, taxonomy_orig))
Niccolai_tax_tab <- Niccolai_tax_tab %>% separate(taxonomy_orig, c('kingdom', 'phylum', 'class', 'order','family', 'genus', 'species'), sep=";")
Niccolai_tax_tab <- Niccolai_tax_tab %>% select(-c(species)) #had to include this for separation's sake but will remove now.
#One might argue that count_tabs have taxonomy already included. but, there are capitalization and other dataset-specific differences
#that make this difficult. 

#import metadata
Burns_sample_info_tab <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Burns_2015_sample_metadata.csv")
Hale_sample_info_tab <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Hale2018/Hale_2018_metadata_for_MICOM_samples.csv")
Niccolai_sample_info_tab <- read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Niccolai2020/Niccolai2020_all_metadata.csv")

#clean up and move on
rm(Burns_abs_otu_table, Burns_count_tab, Hale_abs_otu_table, Hale_count_tab, Niccolai_abs_otu_table, Niccolai_count_tab)

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
Hale_sample_info_tab <- column_to_rownames(Hale_sample_info_tab, var = "sampleIDs_from_read_counts")
#Niccolai_sample_info_tab <- na.omit(Niccolai_sample_info_tab)  #don't run this for Niccolai, as the dataset is a bit sparse
Niccolai_sample_info_tab <- column_to_rownames(Niccolai_sample_info_tab, var = "Sample_ID")

#tax_table tab: set rownames as otu_IDs
Burns_tax_tab <- column_to_rownames(Burns_tax_tab, var="OTU_ID")
Hale_tax_tab <- column_to_rownames(Hale_tax_tab, var="OTU_ID")
Niccolai_tax_tab <- column_to_rownames(Niccolai_tax_tab, var="OTU_ID")

#otu_table/count_tab: set rownamses as otu_id
Burns_count_tab2 <- column_to_rownames(Burns_count_tab2, var="OTU_ID")
Burns_count_tab2 <- Burns_count_tab2 %>% select(-c("taxonomy")) #some reason there's an extra column here
Hale_count_tab2 <- column_to_rownames(Hale_count_tab2, var="OTU_ID")
Niccolai_count_tab2 <- column_to_rownames(Niccolai_count_tab2, var="OTU_ID")
Niccolai_count_tab2 <- Niccolai_count_tab2 %>% select(-c("Domain", "Phylum", "Class", "Order", "Family", "Genus")) #some reason there's an extra column here

#make phyloseq objects
ps_Burns <- phyloseq(otu_table(as.matrix(Burns_count_tab2), taxa_are_rows=TRUE), 
                        sample_data(Burns_sample_info_tab), 
                        tax_table(as.matrix(Burns_tax_tab)))
ps_Hale <- phyloseq(otu_table(as.matrix(Hale_count_tab2), taxa_are_rows=TRUE), 
                     sample_data(Hale_sample_info_tab), 
                     tax_table(as.matrix(Hale_tax_tab)))
ps_Niccolai <- phyloseq(otu_table(as.matrix(Niccolai_count_tab2), taxa_are_rows=TRUE), 
                     sample_data(Niccolai_sample_info_tab), 
                     tax_table(as.matrix(Niccolai_tax_tab)))

#############Step 4: graphing and plotting
#start with some exploratory analyses
plot_bar(ps_Burns, fill = "phylum")
#not that interesting

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
        strip.text.x = element_text(size = 10, color = "black", face = "bold"))

HaleAlpha <- plot_richness(ps_Hale, x="normal_adjacent_or_tumor_tissue_specimen", measures=c("Shannon", "Simpson", "Chao1"), 
                           color="normal_adjacent_or_tumor_tissue_specimen", title="Hale 2018 dataset")+
  geom_boxplot(color="black", aes(fill= normal_adjacent_or_tumor_tissue_specimen)) +
  geom_point(color="black") +
  theme_bw() +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
  theme(legend.position='none', axis.title.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("tumor", "normal")), label = "p.signif") +
  theme(strip.background = element_rect(color="black", size=3, linetype="blank"), strip.placement = "outside",
        strip.text.x = element_text(size = 10, color = "black", face = "bold"))

NiccolaiAlpha <- plot_richness(ps_Niccolai, x="Description", measures=c("Shannon", "Simpson", "Chao1"), color="Description",
                            title="Niccolai 2020 dataset")+
  geom_boxplot(color="black", aes(fill= Description)) +
  geom_point(color="black") +
  theme_bw() +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
  theme(legend.position='none', axis.title.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("tumor", "normal")), label = "p.signif") +
  theme(strip.background = element_rect(color="black", size=3, linetype="blank"), strip.placement = "outside",
        strip.text.x = element_text(size = 10, color = "black", face = "bold"))

#combine graphs
alpha_diversity <- ggarrange(BurnsAlpha,HaleAlpha,NiccolaiAlpha, 
          ncol = 2, nrow = 2)

alpha_diversity

###########Beta Diversity
# Transform data to proportions as appropriate for Bray-Curtis distances
ps_Burns.prop <- transform_sample_counts(ps_Burns, function(otu) otu/sum(otu))
ps_Hale.prop <- transform_sample_counts(ps_Hale, function(otu) otu/sum(otu))
ps_Niccolai.prop <- transform_sample_counts(ps_Niccolai, function(otu) otu/sum(otu))

ord.nmds.bray.Burns <- ordinate(ps_Burns.prop, method="NMDS", distance="bray", trymax=250)
ord.nmds.bray.Hale <- ordinate(ps_Hale.prop, method="NMDS", distance="bray", trymax=350)
ord.nmds.bray.Niccolai <- ordinate(ps_Niccolai.prop, method="NMDS", distance="bray", trymax=200)

ord.PCoA.bray.Burns <- ordinate(ps_Burns.prop, method="PCoA", distance="bray") 
ord.PCoA.bray.Hale <- ordinate(ps_Hale.prop, method="PCoA", distance="bray") 
ord.PCoA.bray.Niccolai <- ordinate(ps_Niccolai.prop, method="PCoA", distance="bray")

plot_ordination(ps_Burns.prop, ord.nmds.bray.Burns, color="Description", title="Burns 2015, Bray NMDS by Description")
plot_ordination(ps_Burns.prop, ord.PCoA.bray.Burns, color="Description", title="Burns 2015, Bray PCoA by Description")

plot_ordination(ps_Hale.prop, ord.nmds.bray.Hale, color="normal_adjacent_or_tumor_tissue_specimen", title="Hale 2018, Bray NMDS by Description")
plot_ordination(ps_Hale.prop, ord.PCoA.bray.Hale, color="normal_adjacent_or_tumor_tissue_specimen", title="Hale 2015, Bray PCoA by Description")

plot_ordination(ps_Niccolai.prop, ord.nmds.bray.Niccolai, color="Description", title="Niccolai 2020, Bray NMDS by Description")
plot_ordination(ps_Niccolai.prop, ord.PCoA.bray.Niccolai, color="Description", title="Niccolai 2020, Bray PCoA by Description")

#none of these ordination plots look nice. 

plot_bar(ps_Burns.prop, x="Description", fill="phylum") + geom_bar(stat="identity") 
plot_bar(ps_Burns, x="Description", fill="phylum") + geom_bar(stat="identity") 

###########Deploying a CLR transform
library(microbiome)

ps_Burns_clr <- microbiome::transform(ps_Burns, "clr")
ps_Hale_clr <- microbiome::transform(ps_Hale, "clr")
ps_Niccolai_clr <- microbiome::transform(ps_Niccolai, "clr")

#ake a quick look
phyloseq::otu_table(ps_Burns_clr)[1:5, 1:5]

#make pca
ord.PCoA.bray.Burns.clr <- ordinate(ps_Burns_clr, method="PCoA", distance="bray") 
ord.PCoA.bray.Hale.clr <- ordinate(ps_Hale_clr, method="PCoA", distance="bray") 
ord.PCoA.bray.Niccolai.clr <- ordinate(ps_Niccolai_clr, method="PCoA", distance="bray")

ps_Burns_clr2 <- phyloseq::ordinate(ps_Burns_clr, "RDA")
ps_Hale_clr2 <- phyloseq::ordinate(ps_Hale_clr, "RDA")
ps_Niccolai_clr2 <- phyloseq::ordinate(ps_Niccolai_clr, "RDA")

#phyloseq::plot_scree(ps_Burns_clr2) 

#let's run permanova
#Generate distance matrix
Burns_clr_dist_matrix <- phyloseq::distance(ps_Burns_clr, method = "euclidean") 
Hale_clr_dist_matrix <- phyloseq::distance(ps_Hale_clr, method = "euclidean") 
Niccolai_clr_dist_matrix <- phyloseq::distance(ps_Niccolai_clr, method = "euclidean") 
#ADONIS test
vegan::adonis(Burns_clr_dist_matrix ~ phyloseq::sample_data(ps_Burns_clr)$Description) #hey, p=0.004
vegan::adonis(Hale_clr_dist_matrix ~ phyloseq::sample_data(ps_Hale_clr)$normal_adjacent_or_tumor_tissue_specimen) # p=1
vegan::adonis(Niccolai_clr_dist_matrix ~ phyloseq::sample_data(ps_Niccolai_clr)$Description) #p=0.01

#The test from ADONIS can be confounded by differences in dispersion (or spread)…so we want to check this as well.
Burns_dispr <- vegan::betadisper(Burns_clr_dist_matrix, phyloseq::sample_data(ps_Burns_clr)$Description)
permutest(Burns_dispr) #p=0.605, so there is no difference in dispersion
Hale_dispr <- vegan::betadisper(Hale_clr_dist_matrix, phyloseq::sample_data(ps_Hale_clr)$normal_adjacent_or_tumor_tissue_specimen)
permutest(Hale_dispr) #p=0.282
Niccolai_dispr <- vegan::betadisper(Niccolai_clr_dist_matrix, phyloseq::sample_data(ps_Niccolai_clr)$Description)
permutest(Niccolai_dispr) #p=0.274

#Now, graph
Burns_pca <- phyloseq::plot_ordination(ps_Burns_clr, ps_Burns_clr2, type="samples", color="Description") + 
  geom_point(size = 2) +
  stat_ellipse(aes(group = Description), linetype = 2)+
  theme_bw() +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
  labs(title="Burns 2015 dataset", subtitle= "PERMANOVA p=0.004") +
  theme(plot.title = element_text(size=14, hjust = 0.5), plot.subtitle=element_text(hjust=0.5))

Hale_pca <- phyloseq::plot_ordination(ps_Hale_clr, ps_Hale_clr2, type="samples", color="normal_adjacent_or_tumor_tissue_specimen") + 
  geom_point(size = 2) +
  stat_ellipse(aes(group = normal_adjacent_or_tumor_tissue_specimen), linetype = 2) +
  theme_bw() +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
  labs(title="Hale 2018 dataset", color="Description", subtitle= "PERMANOVA p>0.99") +
  theme(plot.title = element_text(size=14, hjust = 0.5), plot.subtitle=element_text(hjust=0.5))

Niccolai_pca <- phyloseq::plot_ordination(ps_Niccolai_clr, ps_Niccolai_clr2, type="samples", color="Description") + 
  geom_point(size = 2) +
  stat_ellipse(aes(group = Description), linetype = 2) +
  theme_bw() +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
  labs(title="Niccolai 2020 dataset", subtitle= "PERMANOVA p=0.01") +
  theme(plot.title = element_text(size=14, hjust = 0.5), plot.subtitle=element_text(hjust=0.5))

#combine graphs
pcas <- ggarrange(Burns_pca,Hale_pca,Niccolai_pca, 
                  labels=c("B.", "C.", "D."),
                             ncol = 3, nrow = 1)
pcas #add , legend.position='bottom' to themes to get this to look nice

pcas2 <- ggarrange(Burns_pca,Hale_pca,Niccolai_pca, 
                  labels=c("B.", "C.", "D."),
                  ncol = 1, nrow = 3)
pcas2

#################Step 5: DEseq2
library("DESeq2")

Burns_dds = phyloseq_to_deseq2(ps_Burns, ~factor(Patient_Blind_ID) + factor(Description))
Burns_dds = DESeq(Burns_dds, test="Wald", fitType="parametric")
Burns_res <- results(Burns_dds)
Burns_res <- Burns_res[order(Burns_res$padj),] #re-order by p-value
head(Burns_res)
Burns_sigtab = Burns_res[which(Burns_res$padj < 0.05), ]
#import taxonomic info
Burns_sigtab = cbind(as(Burns_sigtab, "data.frame"), as(tax_table(ps_Burns)[rownames(Burns_sigtab), ], "matrix"))
head(Burns_sigtab)
#No fusobacterium, only weird shit. Let's try the other ones.

Hale_dds = phyloseq_to_deseq2(ps_Hale, ~factor(host_subject_id) + factor(normal_adjacent_or_tumor_tissue_specimen))
Hale_dds = DESeq(Hale_dds, test="Wald", fitType="parametric", sfType = "poscounts") #add poscounts to avoid "every gene contains at least one zero" error
Hale_res <- results(Hale_dds)
Hale_res <- Hale_res[order(Hale_res$padj),] #re-order by p-value
head(Hale_res)
Hale_sigtab = Hale_res[which(Hale_res$padj < 0.05), ]
#import taxonomic info
Hale_sigtab = cbind(as(Hale_sigtab, "data.frame"), as(tax_table(ps_Hale)[rownames(Hale_sigtab), ], "matrix"))
head(Hale_sigtab)

Niccolai_dds = phyloseq_to_deseq2(ps_Niccolai, ~factor(Patient_Blind_ID) + factor(Description))
Niccolai_dds = DESeq(Niccolai_dds, test="Wald", fitType="parametric", sfType = "poscounts")
Niccolai_res <- results(Niccolai_dds)
Niccolai_res <- Niccolai_res[order(Niccolai_res$padj),] #re-order by p-value
head(Niccolai_res)
Niccolai_sigtab = Niccolai_res[which(Niccolai_res$padj < 0.05), ]
#import taxonomic info
Niccolai_sigtab = cbind(as(Niccolai_sigtab, "data.frame"), as(tax_table(ps_Niccolai)[rownames(Niccolai_sigtab), ], "matrix"))
head(Niccolai_sigtab)

#yeah, I"m pretty sure this is completely meaningless.
require("ggrepel")
library(viridis)
theme_set(theme_bw())
ggplot(Burns_sigtab, aes(x=order, y=log2FoldChange, color=order)) + geom_point(size=5) + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust=0.5, size=10), axis.text.y= element_text(size=10)) +
  scale_y_continuous(name = "log2FoldChange in taxon abundance\n normal(base) vs tumor", limits=c(-30,30)) +
  scale_x_discrete(name = "Order", guide = guide_axis(n.dodge=2)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  labs(title = expression("Burns 2015, differentially abundant taxa"), size=14) +
  theme(plot.title = element_text(hjust=0.5), legend.position='none') +
  geom_label_repel(aes(label = signif(padj, digits=3)), color="black", box.padding= 0.5, size=3, xlim = c(NA, Inf)) 

ggplot(Hale_sigtab, aes(x=family, y=log2FoldChange, color=family)) + geom_point(size=5) + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust=0.5, size=10), axis.text.y= element_text(size=10)) +
  scale_y_continuous(name = "log2FoldChange in taxon abundance\n normal(base) vs tumor", limits=c(-30,30)) +
  scale_x_discrete(name = "Family", guide = guide_axis(n.dodge=2)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  labs(title = expression("Hale 2018, differentially abundant taxa"), size=14) +
  theme(plot.title = element_text(hjust=0.5), legend.position='none') +
  geom_label_repel(aes(label = signif(padj, digits=3)), color="black", box.padding= 0.5, size=3, xlim = c(NA, Inf)) 

ggplot(Niccolai_sigtab, aes(x=genus, y=log2FoldChange, color=genus)) + geom_point(size=5) + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust=0.5, size=10), axis.text.y= element_text(size=10)) +
  scale_y_continuous(name = "log2FoldChange in taxon abundance\n normal(base) vs tumor", limits=c(-30,30)) +
  scale_x_discrete(name = "Genus", guide = guide_axis(n.dodge=2)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  labs(title = expression("Niccolai 2020, differentially abundant taxa"), size=14) +
  theme(plot.title = element_text(hjust=0.5), legend.position='none') +
  geom_label_repel(aes(label = signif(padj, digits=3)), color="black", box.padding= 0.5, size=3, xlim = c(NA, Inf)) 

#well.... Burns has too few interesting taxa, and Niccolai has too many.

#can we just look at fuso?
ggplot(subset(phylum="Fusobacteria" %in% Burns_sigtab$phylum), aes(x=order, y=log2FoldChange, color=order)) + geom_point(size=5) + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust=0.5, size=10), axis.text.y= element_text(size=10)) +
  scale_y_continuous(name = "log2FoldChange in taxon abundance\n normal(base) vs tumor", limits=c(-30,30)) +
  scale_x_discrete(name = "Order", guide = guide_axis(n.dodge=2)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  labs(title = expression("Burns 2015, differentially abundant taxa"), size=14) +
  theme(plot.title = element_text(hjust=0.5), legend.position='none') +
  geom_label_repel(aes(label = signif(padj, digits=3)), color="black", box.padding= 0.5, size=3, xlim = c(NA, Inf)) 

ggplot(subset(Hale_sigtab, phylum %in% " Fusobacteria"), aes(x=genus, y=log2FoldChange, color=genus))+ geom_point(size=5)+ 
  theme(axis.text.x = element_text(hjust = 0.5, vjust=0.5, size=10), axis.text.y= element_text(size=10)) +
  scale_y_continuous(name = "log2FoldChange in taxon abundance\n normal(base) vs tumor", limits=c(0,30)) +
  scale_x_discrete(name = "Genus") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  labs(title = expression("Hale 2018, Phylum Fusobacteria"), size=14) +
  theme(plot.title = element_text(hjust=0.5), legend.position='none') +
  geom_label_repel(aes(label = signif(padj, digits=3)), color="black", box.padding= 0.5, size=3, xlim = c(NA, Inf)) 

ggplot(subset(Niccolai_sigtab, phylum %in% " Fusobacteria"), aes(x=genus, y=log2FoldChange, color=genus))+ geom_point(size=5)+ 
  theme(axis.text.x = element_text(hjust = 0.5, vjust=0.5, size=10), axis.text.y= element_text(size=10)) +
  scale_y_continuous(name = "log2FoldChange in taxon abundance\n normal(base) vs tumor", limits=c(0,15)) +
  scale_x_discrete(name = "Genus") +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  labs(title = expression("Niccolai 2020, Phylum Fusobacteria"), size=14) +
  theme(plot.title = element_text(hjust=0.5), legend.position='none') +
  geom_label_repel(aes(label = signif(padj, digits=3)), color="black", box.padding= 0.5, size=3, xlim = c(NA, Inf)) 




