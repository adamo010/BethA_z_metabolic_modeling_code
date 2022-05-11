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
library("DESeq2")
library("BiocManager")
library("microbiome")
library("readxl")
library("forcats")
library("viridis")

#basing this on 12.09.20_MICOM_cooop_tradeoff_growth_rates.R
#further based off 01.25.21_MICOM_coop_tradeoff_EUstd_diet.R; adapted for Hale2018 data from Burns data

#now, I would like to plot dot plots of all the comm GRs, with tradeoff on the x-axis, GR on the y-axis, color by cancer status.
##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Hale2018/") 

##########step 1: import data and clean up
#this dataframe contains the longform data; see 08.31.21_Hale2018_comm_GR_combining_files_and_adding_metadata.py for generating it. 
group_GR_all <- read_csv("08.31.21_Hale2018_V12MSIrun_community_GRs_combined_plus_ranks.csv")

#set theme to have nice colours
theme_set(theme_bw())

#set host_subject_id as factor
group_GR_all$host_subject_id <- as.factor(group_GR_all$host_subject_id)

#rename a few variables
group_GR_all <- group_GR_all %>% rename(Description = normal_adjacent_or_tumor_tissue_specimen)

##########step 2: graphing!
#a terrible graph.
ggplot(group_GR_all, aes(x=host_subject_id, y=comm_gr, colour= Description)) +
  geom_dotplot(aes(fill= Description), binaxis = "y", stackdir= "center", dotsize = 0.4) + 
  ggtitle("Average community growth rates by host_subject_id,\nHale2018 dataset") +
  theme(plot.title = element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average community growth\nrate, mmol/(gDW * h)", limits=c(0,1)) +
  scale_x_discrete(guide=guide_axis(n.dodge=2)) +
  rotate_x_text(angle = 90) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1")

#I don't like this; it is a mess. 
#The above comment is old; reader, this plot is still a mess. 

#now look at paired samples. Also still a mess.
ggplot(group_GR_all, aes(x=Description, y=comm_gr)) +
  geom_boxplot(aes(color= Description)) +
  geom_point(aes(color= Description)) + 
  geom_line(aes(group= host_subject_id)) +
  ggtitle("Average community growth rates by tumor/normal status,\nHale2018 dataset") +
  theme(plot.title = element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average community growth\nrate, mmol/(gDW * h)", limits=c(0,0.8)) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  stat_compare_means(paired= TRUE)

##########step 3: graphing and stats!
#this code sucks. I'm borrowing some code and graphs from 02.17.21_graphing_taxon_specific_GRs_with_stats.R
#since I don't want to make more dataframes in python, let's try working with the same code here.
#BUT I do need to add a "difference_direction" column. So i want to convert to wide form- create a comm_gr_tumor and comm_gr_normal column
group_GR_all_wide = group_GR_all %>% spread(Description, comm_gr) #split comm_gr column into 2 columns, named by different values of Description
group_GR_all_wide <- group_GR_all_wide %>% rename(comm_gr_normal = normal, comm_gr_tumor=tumor) #rename these columns
group_GR_all_wide <- select(group_GR_all_wide, -c(...1, Sample_ID)) #drop the columns called ...1 (the index column, and the sample ID column)
#now, collapse the dataframe by host_subject_id, so there are no blanks in comm_gr_normal and comm_gr_tumor columns
merge_fun <- function(x) x[!is.na(x)]
group_GR_all_wide2 <-group_GR_all_wide %>%
  group_by(host_subject_id) %>%
  summarise_all(funs(merge_fun))
#that, weirdly, duplicated a bunch of rows. Remove these.
group_GR_all_wide2 <- group_GR_all_wide2[!duplicated(group_GR_all_wide2$host_subject_id), ]
#now, create a new column called mean_GR_difference
group_GR_all_wide2 <- group_GR_all_wide2 %>%
  mutate(mean_GR_difference=comm_gr_tumor-comm_gr_normal)
#now, create another new column called Difference_direction
group_GR_all_wide2 <- group_GR_all_wide2 %>%
  mutate(Difference_direction = case_when(
    mean_GR_difference > 0 ~ "tumor > normal",
    mean_GR_difference < 0 ~ "tumor < normal",
    mean_GR_difference == 0 ~ "Zero change",
  ))
#cool. Now let's make some graphs and run some stats.
# Compute t-test- REMEMBER to use the non wide form of the dataframe
comm_gr_by_description <- wilcox.test(comm_gr ~ Description, data = group_GR_all, paired = TRUE)
comm_gr_by_description$p.value #print p-value

#GRAPH!!!
ggpaired(group_GR_all_wide2, cond1= "comm_gr_normal", cond2= "comm_gr_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Community growth rates, Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.990")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average taxon growth\nrate, mmol/(gDW * h)", limits=c(0,0.8)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill='none')



