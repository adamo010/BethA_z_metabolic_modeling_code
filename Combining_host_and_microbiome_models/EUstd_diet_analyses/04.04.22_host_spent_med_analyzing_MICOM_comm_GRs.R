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
library("readxl")
library("forcats")
library("viridis")

#the goal here is to analyze the MICOM community growth rates from communities grown on host spent medium
#start with 11.08.21_Niccolai2020_comm_GR_analysis.R and 03.21.22_combining_host_and_microbiomes_EUstd_host_GR_analysis.R

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_analyses/") 

##########step 1: import data and clean up
#this dataframe contains the longform data; see 08.31.21_Hale2018_comm_GR_combining_files_and_adding_metadata.py for generating it. 
group_GR_all <- read_csv("04.01.22_combining_host_and_microbiomes_EUstd_microbiome_comm_GRs.csv")
group_GR_all <- select(group_GR_all, -c(...1)) #drop the columns called ...1 (the index column, and the sample ID column)
#set theme to have nice colours
theme_set(theme_bw())
#set Patient_Blind_ID as factor
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
  labs(title = expression("MICOM community growth rates, Burns2015 data"), 
       subtitle= expression("p-values from uncorrected pairwise Mann-Whitney U tests")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Community growth\nrate, mmol/(gDW * h)", limits=c(0,0.4)) +
  scale_x_discrete(labels = c('Normal\nmicrobiome,\nnormal host','Normal\nmicrobiome,\ntumor host',
                              'Tumor\nmicrobiome,\nnormal host','Tumor\nmicrobiome,\ntumor host'),
                   name= "Patient microbiome and host tissue descriptions") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "none") +
  guides(fill='none') +
  stat_compare_means(comparisons= y_comparisons, method = "wilcox.test", size= 3, label.y = c(0.20,0.24,0.28,0.32,0.36,0.40))

