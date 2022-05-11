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

#the purpose of this this code is to graph the Fuso-specific exchange fluxes
#borrowed from 11.08.21_Burns2015/11.09.21_Hale2018/11.10.21_Niccolai20202_Fuso_fluxpath_analysis_exchanges_only.R
#from 11.10 version- new fusobacterium filtering pipeline

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/") 

##########import data
Burns_data <- read.csv("Analyzing_MICOM_communities_Burns2015/12.20.21_Fusobacterium_exchange_fluxes_only_alldata_Burns2015_data.csv")
Burns_stats <- read.csv("Analyzing_MICOM_communities_Burns2015/12.20.21_Fusobacterium_exchange_fluxes_only_stats_Burns2015_data.csv")
Hale_data <- read.csv("Analyzing_MICOM_communities_Hale2018/12.21.21_Fusobacterium_exchange_fluxes_only_alldata_Hale2018_data.csv")
Hale_stats <- read.csv("Analyzing_MICOM_communities_Hale2018/12.21.21_Fusobacterium_exchange_fluxes_only_stats_Hale2018_data.csv")
Niccolai_data <- read.csv("Analyzing_MICOM_communities_Niccolai2020/12.21.21_Fusobacterium_exchange_fluxes_only_alldata_Niccolai2020_data.csv")
Niccolai_stats <- read.csv("Analyzing_MICOM_communities_Niccolai2020/12.21.21_Fusobacterium_exchange_fluxes_only_stats_Niccolai2020_data.csv")

#let's see if any of the top... 10 metabolites, by q-value, are in common
Burns_stats <- Burns_stats[order(Burns_stats$by_flux_qvals),] #re-order by q-value
Hale_stats <- Hale_stats[order(Hale_stats$by_flux_qvals),] #re-order by q-value
Niccolai_stats <- Niccolai_stats[order(Niccolai_stats$by_flux_qvals),] #re-order by q-value

Burns_stats2 <- Burns_stats %>% top_n(10)
Hale_stats2 <- Hale_stats %>% top_n(10)
Niccolai_stats2 <- Niccolai_stats %>% top_n(10)

#thiamin is in all three, as is uracil and water. None of them are even slightly significant. 


########Graphing proline
#Proline, Burns
ggpaired(subset(Burns_data, fluxpath_name %in% c("EX_pro_L(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = "Fusobacterium proline exchange flux,\n Burns2015 data", 
       subtitle = expression("paired Wilcoxan, p=0.0452, q= 0.217")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Sum of Fusobacterium fluxes,\nmmol/(gDW * h)", limits=c(-200,200)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Proline, Hale
ggpaired(subset(Hale_data, fluxpath_name %in% c("EX_pro_L(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = "Fusobacterium proline exchange flux,\n Hale2018 data", 
       subtitle = expression("paired Wilcoxan, p=0.0164, q= 0.0982")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Sum of Fusobacterium fluxes,\nmmol/(gDW * h)", limits=c(-700,300)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Proline, Niccolai
ggpaired(subset(Niccolai_data, fluxpath_name %in% c("EX_pro_L(e)")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = "Fusobacterium proline exchange flux,\n Niccolai2020 data", 
       subtitle = expression("paired Wilcoxan, p=0.0301, q= 0.166")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5), plot.margin = unit(c(1,0.1,0.1,0.1), "cm")) +
  scale_y_continuous(name = "Sum of Fusobacterium fluxes,\nmmol/(gDW * h)", limits=c(-500,200)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")




