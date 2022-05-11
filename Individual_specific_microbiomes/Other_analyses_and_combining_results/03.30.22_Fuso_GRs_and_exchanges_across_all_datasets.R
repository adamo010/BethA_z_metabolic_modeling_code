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
library("ggtext")

#The goal of this script is to centralize the Fusobacterium data from Burns, Hale, and Niccolai datasets. We want to graph growth rates
#and exchange fluxes.
#the goal isn't to re-analyze; just to graph.

#the growth rate source data for this file is:
#Burns: 09.24.21_Burns2015_indiv_spp_GR_analysis.R; 10.27.21_identifying_fuso_containing_Burns_samples.R
#Hale: 09.23.21_Hale2018_indiv_spp_GR_analysis.R
#Niccolai: 11.08.21_Niccolai2020_indiv_spp_GR_analysis

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Other_analyses_and_combining_results/") 

##########Importing data and extracting Fuso only
Burns_gr_data <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/09.24.21_filtered_microbial_GR_data_all_Burns2015_data.csv")
Burns_gr_data <- subset(Burns_gr_data, select=-c(X,...1))
Burns_fuso_gr_data <- Burns_gr_data[(Burns_gr_data$genus== "Fusobacterium"),]
Burns_gr_stats <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/09.24.21_filtered_microbial_GR_all_results_Burns2015_data.csv")
Burns_gr_stats <- subset(Burns_gr_stats, select=-c(X))
Burns_fuso_gr_stats <- Burns_gr_stats[(Burns_gr_stats$genus_ids_point_one=="Fusobacteriaceae_Fusobacterium"),]
rm(Burns_gr_data, Burns_gr_stats)

Hale_gr_data <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Hale2018/09.23.21_filtered_microbial_GR_data_all_Hale2018_data.csv")
Hale_gr_data <- subset(Hale_gr_data, select=-c(X, ...1))
Hale_fuso_gr_data <- Hale_gr_data[(Hale_gr_data$Genus== "Fusobacterium"),]
Hale_gr_stats <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Hale2018/09.23.21_filtered_microbial_GR_all_results_Hale2018_data.csv")
Hale_gr_stats <- subset(Hale_gr_stats, select=-c(X))
Hale_fuso_gr_stats <- Hale_gr_stats[(Hale_gr_stats$genus_ids_point_one=="Fusobacteriaceae_Fusobacterium"),]
rm(Hale_gr_data, Hale_gr_stats)

Niccolai_gr_data <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Niccolai2020/11.08.21_filtered_microbial_GR_data_all_Niccolai2020_data.csv")
Niccolai_gr_data <- subset(Niccolai_gr_data, select=-c(X, ...1))
Niccolai_fuso_gr_data <- Niccolai_gr_data[(Niccolai_gr_data$Genus== "Fusobacterium"),]
Niccolai_gr_stats <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Niccolai2020/11.08.21_filtered_microbial_GR_all_results_Niccolai2020_data.csv")
Niccolai_gr_stats <- subset(Niccolai_gr_stats, select=-c(X))
Niccolai_fuso_gr_stats <- Niccolai_gr_stats[(Niccolai_gr_stats$genus_ids_point_one=="Fusobacteriaceae_Fusobacterium"),]
rm(Niccolai_gr_data, Niccolai_gr_stats)

##########Making graphs
#Burns
Burns_gr_graph <- ggpaired(Burns_fuso_gr_data, cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Burns2015 data"), subtitle = expression("p=1.6x10"^-3*", q= 0.041")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "*Fusobacterium* growth  \nrate, mmol/(gDW * h)", limits=c(0,0.3)) + #need two spaces for line break in markdown
  scale_x_discrete(name= "", labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2", labels = c("tumor < normal", "tumor > normal", "tumor = normal")) +
  theme(legend.position = "none", legend.title=element_blank(), legend.text=element_text(size=11), 
        axis.title.y = ggtext::element_markdown(), legend.key.width= unit(1.5, 'cm')) +
  guides(fill="none", color=guide_legend(nrow=1, byrow=TRUE, override.aes = list(size=2)))

#Hale
Hale_gr_graph <- ggpaired(Hale_fuso_gr_data, cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Hale2018 data"), subtitle = expression("p=3.6x10"^-4*", q= 0.011")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "", limits=c(0,0.6)) +
  scale_x_discrete(name= "", labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2", labels = c("tumor < normal", "tumor > normal", "tumor = normal")) +
  theme(legend.position = "none", legend.title=element_blank(), legend.text=element_text(size=11), 
        axis.title.y = ggtext::element_markdown(), legend.key.width= unit(1.5, 'cm')) +
  guides(fill="none", color=guide_legend(nrow=1, byrow=TRUE, override.aes = list(size=2)))

#Niccolai
Niccolai_gr_graph <- ggpaired(Niccolai_fuso_gr_data, cond1= "mean_genus_GR_normal", cond2= "mean_genus_GR_tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Niccolai2020 data"), subtitle = expression("p=2.3x10"^-6*", q= 8.3x10"^-5)) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "", limits=c(0,0.6)) +
  scale_x_discrete(name= "", labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2", labels = c("tumor < normal", "tumor > normal", "tumor = normal")) +
  theme(legend.position = "none", legend.title=element_blank(), legend.text=element_text(size=11), 
        axis.title.y = ggtext::element_markdown(), legend.key.width= unit(1.5, 'cm')) +
  guides(fill="none", color=guide_legend(nrow=1, byrow=TRUE, override.aes = list(size=2)))

#great. Now we have all our graphs. Let's make them into a figure using ggarrange: a NEW command
all_grs_fig <- ggarrange(Burns_gr_graph, Hale_gr_graph, Niccolai_gr_graph,
                    labels = c("D.", "E.", "F."),
                    common.legend = TRUE, legend = "bottom",
                    ncol = 3, nrow = 1)
all_grs_fig #width by height 1000 X 500
#this is an extremely useful piece of code

###################################################################OKAY! Next: fluxes.
#for all samples, we're going with 
#For fuso only, we're going with Valine exchange EX_val_L(e) because it's the only metabolite that showed up in the top 20 for all three datasets and is *kind of* significant.
#for all datasets, co2 exchange is the only shared metab, and it's not sig. 
#the growth rate source data for this file is:
#Burns: 10.21.21_BUrns2015_fluxpath_analysis_exchanges_only.R
#Hale: 09.23.21_Hale2018_indiv_spp_GR_analysis.R
#Niccolai: 11.08.21_Niccolai2020_indiv_spp_GR_analysis

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Other_analyses_and_combining_results/") 

##########Importing data: FUSO SPECIFIC FLUXES
Burns_flux_data <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/12.20.21_Fusobacterium_exchange_fluxes_only_alldata_Burns2015_data.csv")
Burns_flux_data <- subset(Burns_flux_data, select=-c(X))
Burns_val_flux_data <- Burns_flux_data[(Burns_flux_data$fluxpath_name== "EX_val_L(e)"),]
Burns_flux_stats <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/12.20.21_Fusobacterium_exchange_fluxes_only_stats_Burns2015_data.csv")
Burns_flux_stats <- subset(Burns_flux_stats, select=-c(X))
Burns_filtered_flux_stats <- Burns_flux_stats[(Burns_flux_stats$fluxpath_name=="EX_val_L(e)"),]
rm(Burns_flux_data, Burns_flux_stats)

Hale_flux_data <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Hale2018/12.21.21_Fusobacterium_exchange_fluxes_only_alldata_Hale2018_data.csv")
Hale_flux_data <- subset(Hale_flux_data, select=-c(X))
Hale_val_flux_data <- Hale_flux_data[(Hale_flux_data$fluxpath_name== "EX_val_L(e)"),]
Hale_flux_stats <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Hale2018/12.21.21_Fusobacterium_exchange_fluxes_only_stats_Hale2018_data.csv")
Hale_flux_stats <- subset(Hale_flux_stats, select=-c(X))
Hale_fuso_flux_stats <- Hale_flux_stats[(Hale_flux_stats$fluxpath_name=="EX_val_L(e)"),]
rm(Hale_flux_data, Hale_flux_stats)

Niccolai_flux_data <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Niccolai2020/12.21.21_Fusobacterium_exchange_fluxes_only_alldata_Niccolai2020_data.csv")
Niccolai_flux_data <- subset(Niccolai_flux_data, select=-c(X))
Niccolai_val_flux_data <- Niccolai_flux_data[(Niccolai_flux_data$fluxpath_name== "EX_val_L(e)"),]
Niccolai_flux_stats <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Niccolai2020/12.21.21_Fusobacterium_exchange_fluxes_only_stats_Niccolai2020_data.csv")
Niccolai_flux_stats <- subset(Niccolai_flux_stats, select=-c(X))
Niccolai_fuso_flux_stats <- Niccolai_flux_stats[(Niccolai_flux_stats$fluxpath_name=="EX_val_L(e)"),]
rm(Niccolai_flux_data, Niccolai_flux_stats)

##########Making graphs
#Burns
Burns_flux_graph <- ggpaired(Burns_val_flux_data, cond1= "normal", cond2= "tumor",
                           y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Fusobacterium valine exchange, Burns2015 data"), subtitle = expression("paired Wilcoxan, p=0.010, q= 0.0626")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Exchange flux,\nmmol/(gDW * h)", limits=c(-500,500)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Hale
Hale_flux_graph <- ggpaired(Hale_val_flux_data, cond1= "normal", cond2= "tumor",
                             y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Fusobacterium valine exchange, Hale2018 data"), subtitle = expression("paired Wilcoxan, p=0.0752, q= 0.301")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Exchange flux,\nmmol/(gDW * h)", limits=c(-10,200)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")


#Niccolai
Niccolai_flux_graph <- ggpaired(Niccolai_val_flux_data, cond1= "normal", cond2= "tumor",
                            y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Fusobacterium valine exchange, Niccolai2020 data"), subtitle = expression("paired Wilcoxan, p=0.00182, q= 0.0697")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Exchange flux,\nmmol/(gDW * h)", limits=c(-200,600)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#great. Now we have all our graphs. Let's make them into a figure using ggarrange: a NEW command
all_flux_fig <- ggarrange(Burns_flux_graph, Hale_flux_graph, Niccolai_flux_graph,
                         labels = c("D.", "E.", "F."),
                         ncol = 3, nrow = 1)

###############################I would also like to make some venn diagrams.
if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
library(ggVennDiagram)

#re-import all/unfiltered flux stats and filter by top 20
Burns_flux_stats <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/12.20.21_Fusobacterium_exchange_fluxes_only_stats_Burns2015_data.csv")
Burns_flux_stats <- subset(Burns_flux_stats, select=-c(X)) #delete extra column
Burns_flux_stats_top20 <- Burns_flux_stats %>% top_n(-20, by_flux_pvals) #keep smallest 20 pvals

Hale_flux_stats <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Hale2018/12.21.21_Fusobacterium_exchange_fluxes_only_stats_Hale2018_data.csv")
Hale_flux_stats <- subset(Hale_flux_stats, select=-c(X))
Hale_flux_stats_top20 <- Hale_flux_stats %>% top_n(-20, by_flux_pvals)

Niccolai_flux_stats <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Niccolai2020/12.21.21_Fusobacterium_exchange_fluxes_only_stats_Niccolai2020_data.csv")
Niccolai_flux_stats <- subset(Niccolai_flux_stats, select=-c(X))
Niccolai_flux_stats_top20 <- Niccolai_flux_stats %>% top_n(-20, by_flux_pvals)

#need to merge all datasets
#first, add a column to each dataset that represents "dataset"
Burns_flux_stats_top20 <- Burns_flux_stats_top20 %>% mutate(dataset= "Burns")
Hale_flux_stats_top20 <- Hale_flux_stats_top20 %>% mutate(dataset= "Hale")
Niccolai_flux_stats_top20 <- Niccolai_flux_stats_top20 %>% mutate(dataset= "Niccolai")
#then, stack dataframes
merged_flux_stats <-rbind(Burns_flux_stats_top20, Hale_flux_stats_top20, Niccolai_flux_stats_top20)

#drop unnecessary columns
merged_flux_stats_minimal <- subset(merged_flux_stats, select=-c(fluxpath_name, by_flux_pvals, by_flux_qvals))

#I think what I need to do is convert this to counts.
#categories of interest: Burns, Hale, Niccolai, Burns-Hale, Burns-Niccolai, Hale-Niccolai, Burns-Hale-Niccolai
merged_flux_stats_counts <- aggregate(dataset ~ fluxpath_description, data= merged_flux_stats_minimal, paste, collapse="_" )
#holy shit, nice. 
#FUCK IT. Venn diagrams... are impossibnle oasdbujifg npaousenjbfpoAUIRWEHBNSGFPOAJSDFBNUGPOVJABNESRDF

merged_flux_stats_counts_t <- transpose(merged_flux_stats_counts)
#redefine row and column names
rownames(merged_flux_stats_counts_t) <- colnames(merged_flux_stats_counts)
colnames(merged_flux_stats_counts_t) <- rownames(merged_flux_stats_counts)

flux_list=as.list(merged_flux_stats_counts_t)
  
merged_flux_stats_counts2 <- table(merged_flux_stats_counts$dataset)

  amerged_flux_stats_counts %>%  #create a new dataframe by filtering the old one
  count(dataset) %>% #
  summarize(n = sum(count_of_nonzeros)) %>% #sums the count_of_nonzeros values within each fluxpath_subsystem and saves the value as n
  mutate(n) #add n as a new column in the dataframe
#create a new dataframe that merges the old and new dataframes





#Then, convert to a list of lists for graphing:
fluxpath_list <- list(
  merged_flux_stats_minimal %>% filter(dataset== "Burns") %>% select(fluxpath_description), 
  merged_flux_stats_minimal %>% filter(dataset== "Hale") %>% select(fluxpath_description), 
  merged_flux_stats_minimal %>% filter(dataset== "Niccolai") %>% select(fluxpath_description),)

#graph:
ggVennDiagram(x <- list(
  merged_flux_stats %>% filter(dataset= "Burns") %>% select(fluxpath_description) %>% unlist() , 
  merged_flux_stats %>% filter(adataset= "Hale") %>% select(fluxpath_description) %>% unlist() , 
  merged_flux_stats %>% filter(dataset= "Niccolai") %>% select(fluxpath_description) %>% unlist(),))

########test data for troubleshooting:
set.seed(20190708)
genes <- paste("gene",1:1000,sep="")
x <- list(
  A = sample(genes,300), 
  B = sample(genes,525), 
  C = sample(genes,440),
  D = sample(genes,350)
)
