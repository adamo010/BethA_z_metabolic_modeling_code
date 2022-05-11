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

#the goal is to graph uptakes of sample-specific CORDA models GROWN ON MICOM SPENT MED to see... differences
#there are three datasets to look at; uptake, secretion, and fluxes. This file is for fluxes.

#generated these data in 03.21.22_combining_MICOM_and_CORDA_output_from_spent_med_V1.py

#starting point is 02.22.22_Burns2022_CORDA_model_analyzing_EUstd_GRs.R

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_analyses/") 

##########Step 2: import data and adjust factors
group_GR_all <- read.csv(file= "03.21.22_combining_host_and_microbiomes_EUstd_host_GRs.csv", header=TRUE)
group_GR_all <- subset(group_GR_all, select=-c(X))
group_GR_all$Patient_Blind_ID <- as.factor(group_GR_all$Patient_Blind_ID)

##########Step 3: stats? I guess?
#what we want to do here is compare GRs within combo_types. So, multiple correct by... six. Six comparisons. Ah ha ha.

#hey, this is fun- test for normality
#Test each group for normality
group_GR_all %>%
  group_by(combo_type) %>%
  summarise(`W Stat` = shapiro.test(growth_rate)$statistic,
            p.value = shapiro.test(growth_rate)$p.value)
#so, the p-values are all below 0.05, which means we should reject the assumption of normality and use the mann-whitney U test

#need to subset data by pairs of types. I REALLY don't want to create new dataframes, but I will...
#UGH FINE FUCK IT.
GRs_set1 <- group_GR_all[ which(group_GR_all$combo_type=='tumor_tumor' | group_GR_all$combo_type== "tumor_normal"), ]
GRs_set2 <- group_GR_all[ which(group_GR_all$combo_type=='tumor_tumor' | group_GR_all$combo_type== "normal_tumor"), ]
GRs_set3 <- group_GR_all[ which(group_GR_all$combo_type=='tumor_tumor' | group_GR_all$combo_type== "normal_normal"), ]
GRs_set4 <- group_GR_all[ which(group_GR_all$combo_type=='tumor_normal' | group_GR_all$combo_type== "normal_tumor"), ]
GRs_set5 <- group_GR_all[ which(group_GR_all$combo_type=='tumor_normal' | group_GR_all$combo_type== "normal_normal"), ]
GRs_set6 <- group_GR_all[ which(group_GR_all$combo_type=='normal_tumor' | group_GR_all$combo_type== "normal_normal"), ]
#the pipe here means 'or'

GR_by_description_set1 <- wilcox.test(growth_rate ~ combo_type, data = GRs_set1, paired = TRUE, exact=FALSE, conf.int=TRUE)
set1_pvalue = GR_by_description_set1$p.value
GR_by_description_set2 <- wilcox.test(growth_rate ~ combo_type, data = GRs_set2, paired = TRUE, exact=FALSE, conf.int=TRUE)
set2_pvalue = GR_by_description_set1$p.value
GR_by_description_set3 <- wilcox.test(growth_rate ~ combo_type, data = GRs_set3, paired = TRUE, exact=FALSE, conf.int=TRUE)
set3_pvalue = GR_by_description_set1$p.value
GR_by_description_set4 <- wilcox.test(growth_rate ~ combo_type, data = GRs_set4, paired = TRUE, exact=FALSE, conf.int=TRUE)
set4_pvalue = GR_by_description_set1$p.value
GR_by_description_set5 <- wilcox.test(growth_rate ~ combo_type, data = GRs_set5, paired = TRUE, exact=FALSE, conf.int=TRUE)
set5_pvalue = GR_by_description_set1$p.value
GR_by_description_set6 <- wilcox.test(growth_rate ~ combo_type, data = GRs_set6, paired = TRUE, exact=FALSE, conf.int=TRUE)
set6_pvalue = GR_by_description_set1$p.value

#ALL. I. DENTICAL!!!
#what if we did it by t-test? also the same

#well, let's graph.
#here are the comparisons we want to do
y_comparisons <- list( c("tumor_tumor", "tumor_normal"), c("tumor_tumor", "normal_tumor"), c("tumor_tumor", "normal_normal"),
                       c("tumor_normal", "normal_tumor"), c("tumor_normal", "normal_normal"), c("normal_tumor", "normal_normal"))
#graph
ggplot(group_GR_all, aes(x= combo_type, y= growth_rate, color=combo_type)) +
  geom_violin(aes(x= combo_type, y= growth_rate, fill=combo_type)) + 
  geom_dotplot(binaxis='y', stackdir='center', stackratio=0.75, dotsize=0.6, aes(x= combo_type, y= growth_rate)) +
  labs(title = expression("Colonocyte growth rates, Burns2015 data"), 
       subtitle= expression("p-values from uncorrected pairwise Mann-Whitney U tests")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Colonocyte growth\nrate, mmol/(gDW * h)", limits=c(0,420)) +
  scale_x_discrete(labels = c('Normal\nmicrobiome,\nnormal host','Normal\nmicrobiome,\ntumor host',
                              'Tumor\nmicrobiome,\nnormal host','Tumor\nmicrobiome,\ntumor host'),
                   name= "Patient microbiome and host tissue descriptions") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "none") +
  guides(fill='none') +
  stat_compare_means(comparisons= y_comparisons, method = "wilcox.test", size= 3, label.y = c(260,290,320,350,380,410))



#################OLD
#import metadata
metadata= read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/CRC_metadata_rnaseq_edited_by_BA.csv")
metadata_slimmed <- select(metadata, c(Tissue_Tube_ID, Patient_Blind_ID, SampleID, Description)) #drop the columns called ...1 (the index column, and the sample ID column)

#merge metadata with data
analyzed_GRs_plus_metadata <-left_join(group_GR_all, metadata_slimmed, by = c("sample_id" = "Tissue_Tube_ID"))

#set factors
analyzed_GRs_plus_metadata$Patient_Blind_ID <- as.factor(analyzed_GRs_plus_metadata$Patient_Blind_ID)
analyzed_GRs_plus_metadata <- select(analyzed_GRs_plus_metadata, -c(SampleID, sample_id))

###########Step 3: rearranging
#split GR column into 2 columns, named by different values of Description
analysed_GRs_wide = analyzed_GRs_plus_metadata %>% spread(Description, growth_rate) 
#now, create a new column called mean_GR_difference
analysed_GRs_wide <- analysed_GRs_wide %>%
  mutate(GR_difference=tumor-normal)
#now, create another new column called Difference_direction
analysed_GRs_wide <- analysed_GRs_wide %>%
  mutate(Difference_direction = case_when(
    GR_difference > 0 ~ "tumor > normal",
    GR_difference < 0 ~ "tumor < normal",
    GR_difference == 0 ~ "Zero change",
  ))

############Step 4: stats and graphing
#cool. Now let's make some graphs and run some stats.
# Compute t-test- REMEMBER to use the non wide form of the dataframe
#GR_by_description <- wilcox.test(growth_rate ~ Description, data = analyzed_GRs_plus_metadata, paired = TRUE)
#can't use wilcoxan- there are pairs with differences=0, so it won't work. 
GR_by_description <- t.test(growth_rate ~ Description, data = analyzed_GRs_plus_metadata, paired = TRUE)
GR_by_description$p.value #print p-value #aha- 0.0205

############Step 5: graph
ggpaired(analysed_GRs_wide, cond1= "normal", cond2= "tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Colonocyte growth rates, Burns2015 data"), subtitle = expression("paired T-test, p=0.0205")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Colonocyte growth\nrate, mmol/(gDW * h)", limits=c(0,300)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill='none')





