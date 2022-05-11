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
library("ggVennDiagram")
library("VennDiagram")

#the goal of this script is to create venn diagrams showing the overlap of tumor vs normal taxa(growth rates) and exchange pathways (fluxes)
#across the Burns/Hale/Niccolai datasets.

#borrow some code from 03.30.22_Fuso_GRs_and_exchanges_across_datasets.R

#clearning out the plots window: dev.off(dev.list()["RStudioGD"]) 

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Other_analyses_and_combining_results/") 

##########Importing data 
#NOTE here that we've already filtered by pathways active in at least half of samples. May need to edit later to filter by significant data.

Burns_gr_data <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/09.24.21_filtered_microbial_GR_data_all_Burns2015_data.csv")
Burns_gr_data <- subset(Burns_gr_data, select=-c(X,...1))
Burns_gr_stats <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/09.24.21_filtered_microbial_GR_all_results_Burns2015_data.csv")
Burns_gr_stats <- subset(Burns_gr_stats, select=-c(X))

Hale_gr_data <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Hale2018/09.23.21_filtered_microbial_GR_data_all_Hale2018_data.csv")
Hale_gr_data <- subset(Hale_gr_data, select=-c(X, ...1))
Hale_gr_stats <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Hale2018/09.23.21_filtered_microbial_GR_all_results_Hale2018_data.csv")
Hale_gr_stats <- subset(Hale_gr_stats, select=-c(X))

Niccolai_gr_data <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Niccolai2020/11.08.21_filtered_microbial_GR_data_all_Niccolai2020_data.csv")
Niccolai_gr_data <- subset(Niccolai_gr_data, select=-c(X, ...1))
Niccolai_gr_stats <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Niccolai2020/11.08.21_filtered_microbial_GR_all_results_Niccolai2020_data.csv")
Niccolai_gr_stats <- subset(Niccolai_gr_stats, select=-c(X))

############Making a venn diagram
#Burns <- unique(Burns_gr_data$genus)
#Hale <- unique(Hale_gr_data$Genus) 
#Niccolai <- unique(Niccolai_gr_data$Genus)

#pick some colours for venn diagram
library("RColorBrewer")
brewer.pal(n = 8, name = "Dark2") # prints "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"

# Helper function to display Venn diagram
#display_venn <- function(x, ...){
  #library(VennDiagram)
  #grid.newpage()
  #venn_object <- venn.diagram(x, filename = NULL, ...)
  #grid.draw(venn_object)
#}

#display_venn(list(Burns=Burns, Hale=Hale, Niccolai=Niccolai),
             #category.names = c("Burns 2015" , "Hale 2018" , "Niccolai 2020"),       
             #fill = c("#1B9E77", "#D95F02", "#7570B3"),
             #alpha = c(0.5, 0.5, 0.5),
             #height = 600, 
             #width = 600,
             #resolution = 300,
             #main= "Number of taxa growing in at least half of patients",
             #main.cex=1.5, cat.cex = 1.2, cex=1.2,
             #cat.default.pos = "outer",
             #cat.pos = c(-27, 27, 0),
             #cat.dist = c(0.055, 0.055, 0.055),
             #main.fontfamily = "sans",
             #cat.fontfamily = "sans",
             #fontfamily = "sans")
#resize to width=550, height=500

#can I get a table of the overlapping taxa?
#threedatasettaxa <- Reduce(intersect, list(Burns, Hale, Niccolai)) #these are the seven taxa that overlap in all three datasets.
#Now, I'd like a column that is a subset of p-value results.
Burns2 <- unique(Burns_gr_data$famgen_col)
Hale2 <- unique(Hale_gr_data$famgen_col) 
Niccolai2 <- unique(Niccolai_gr_data$famgen_col)
threedatasettaxa2 <- Reduce(intersect, list(Burns2, Hale2, Niccolai2)) #these are the seven taxa that overlap in all three datasets.
#now why the hell are the two different versions of this different?
#for fuck's sake. stats are based on famgen, so I guess we should use that.

venn_object_grs <- venn.diagram(list(Burns=Burns2, Hale=Hale2, Niccolai=Niccolai2), filename=NULL, disable_logging=TRUE,
             category.names = c("Burns 2015" , "Hale 2018" , "Niccolai 2020"),       
             fill = c("#1B9E77", "#D95F02", "#7570B3"),
             alpha = c(0.5, 0.5, 0.5),
             height = 600, 
             width = 600,
             resolution = 300,
             main= "No. taxa growing in\n at least half of patients",
             main.cex=1.5, cat.cex = 1.2, cex=1.2,
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 0),
             cat.dist = c(0.055, 0.055, 0.055),
             main.fontfamily = "sans",
             cat.fontfamily = "sans",
             fontfamily = "sans")
grid.draw(venn_object_grs)

#okay. Filter each stats dataframe by what's in threedatasettaxa2
Burns_stats_filt <- Burns_gr_stats[Burns_gr_stats$genus_ids_point_one %in% threedatasettaxa2,]
Burns_stats_filt <- rename(Burns_stats_filt, c(Family_Genus_name=genus_ids_point_one, Burns_p_values=genus_0.1_pvals, Burns_q_values=genus_0.1_qvals))
Burns_stats_filt <- Burns_stats_filt %>% mutate_if(is.numeric, round, 3)
Hale_stats_filt <- Hale_gr_stats[Hale_gr_stats$genus_ids_point_one %in% threedatasettaxa2,]
Hale_stats_filt <- rename(Hale_stats_filt, c(Family_Genus_name=genus_ids_point_one, Hale_p_values=genus_0.1_pvals, Hale_q_values=genus_0.1_qvals))
Hale_stats_filt <- Hale_stats_filt %>% mutate_if(is.numeric, round, 3)
Niccolai_stats_filt <- Niccolai_gr_stats[Niccolai_gr_stats$genus_ids_point_one %in% threedatasettaxa2,]
Niccolai_stats_filt <- rename(Niccolai_stats_filt, c(Family_Genus_name=genus_ids_point_one, Niccolai_p_values=genus_0.1_pvals, Niccolai_q_values=genus_0.1_qvals))
Niccolai_stats_filt <- Niccolai_stats_filt %>% mutate_if(is.numeric, round, 3)

#put all data frames into list
stats_df_list <- list(Burns_stats_filt, Hale_stats_filt, Niccolai_stats_filt)
#merge all data frames in list
stats_df <- stats_df_list %>% reduce(full_join, by='Family_Genus_name')

library(gridExtra)
library(grid)
gr_stats_table <- grid.table(stats_df)

#I don't know how to merge the venn diagram and the table into a single figure. I've spent hours on this. FUCK IT. 
grid.arrange(tableGrob(stats_df),grid.draw(venn_object_grs), nrow=1)

ggarrange(gr_stats_table, grid.draw(venn_object_grs),
                          labels = c("A.","B."),
                          ncol = 2, nrow = 1)


############################## Fluxes
Burns_flux_data <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/12.20.21_Fusobacterium_exchange_fluxes_only_alldata_Burns2015_data.csv")
Burns_flux_data <- subset(Burns_flux_data, select=-c(X))
Burns_flux_stats <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/12.20.21_Fusobacterium_exchange_fluxes_only_stats_Burns2015_data.csv")
Burns_flux_stats <- subset(Burns_flux_stats, select=-c(X)) #delete extra column

Hale_flux_data <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Hale2018/12.21.21_Fusobacterium_exchange_fluxes_only_alldata_Hale2018_data.csv")
Hale_flux_data <- subset(Hale_flux_data, select=-c(X))
Hale_flux_stats <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Hale2018/12.21.21_Fusobacterium_exchange_fluxes_only_stats_Hale2018_data.csv")
Hale_flux_stats <- subset(Hale_flux_stats, select=-c(X))

Niccolai_flux_data <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Niccolai2020/12.21.21_Fusobacterium_exchange_fluxes_only_alldata_Niccolai2020_data.csv")
Niccolai_flux_data <- subset(Niccolai_flux_data, select=-c(X))
Niccolai_flux_stats <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Niccolai2020/12.21.21_Fusobacterium_exchange_fluxes_only_stats_Niccolai2020_data.csv")
Niccolai_flux_stats <- subset(Niccolai_flux_stats, select=-c(X))

#pick some colours for venn diagram
library("RColorBrewer")
brewer.pal(n = 8, name = "Dark2") # prints "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"

#figure out overlaps between data
Burns2 <- unique(Burns_flux_data$fluxpath_name)
Hale2 <- unique(Hale_flux_data$fluxpath_name) 
Niccolai2 <- unique(Niccolai_flux_data$fluxpath_name)
threedatasettaxa2 <- Reduce(intersect, list(Burns2, Hale2, Niccolai2)) #these are the 36 pathways.

venn_object_flux <- venn.diagram(list(Burns=Burns2, Hale=Hale2, Niccolai=Niccolai2), filename=NULL, disable_logging=TRUE,
                            category.names = c("Burns 2015" , "Hale 2018" , "Niccolai 2020"),       
                            fill = c("#1B9E77", "#D95F02", "#7570B3"),
                            alpha = c(0.5, 0.5, 0.5),
                            height = 600, 
                            width = 600,
                            resolution = 300,
                            main= "No. fluxpaths active\n in all patients",
                            main.cex=1.5, cat.cex = 1.2, cex=1.2,
                            cat.default.pos = "outer",
                            cat.pos = c(-27, 27, 7),
                            cat.dist = c(0.1, 0.1, 0.1),
                            main.fontfamily = "sans",
                            cat.fontfamily = "sans",
                            fontfamily = "sans")
grid.draw(venn_object_flux)

dev.off(dev.list()["RStudioGD"]) # Apply dev.off() & dev.list()

#okay. Filter each stats dataframe by what's in threedatasettaxa2
Burns_stats_filt <- Burns_flux_stats[Burns_flux_stats$fluxpath_name %in% threedatasettaxa2,]
Burns_stats_filt <- rename(Burns_stats_filt, c(Burns_p_values=by_flux_pvals, Burns_q_values=by_flux_qvals))
Burns_stats_filt <- Burns_stats_filt %>% mutate_if(is.numeric, round, 3)
Hale_stats_filt <- Hale_flux_stats[Hale_flux_stats$fluxpath_name %in% threedatasettaxa2,]
Hale_stats_filt <- rename(Hale_stats_filt, c(Hale_p_values=by_flux_pvals, Hale_q_values=by_flux_qvals))
Hale_stats_filt <- Hale_stats_filt %>% mutate_if(is.numeric, round, 3)
Niccolai_stats_filt <- Niccolai_flux_stats[Niccolai_flux_stats$fluxpath_name %in% threedatasettaxa2,]
Niccolai_stats_filt <- rename(Niccolai_stats_filt, c(Niccolai_p_values=by_flux_pvals, Niccolai_q_values=by_flux_qvals))
Niccolai_stats_filt <- Niccolai_stats_filt %>% mutate_if(is.numeric, round, 3)

#put all data frames into list
stats_df_list <- list(Burns_stats_filt, Hale_stats_filt, Niccolai_stats_filt)
#merge all data frames in list
stats_df <- stats_df_list %>% reduce(full_join, by='fluxpath_description')
stats_df <- select(stats_df, -c(fluxpath_name.x, fluxpath_name.y, fluxpath_name))

#again, fuck if I know what the right answer is here for making a figure. 
#let's just ignore the tables.

venn_fig <- ggarrange(venn_object_grs, venn_object_flux,
                          labels = c("B.", "C."),
                          ncol = 2, nrow = 1)
venn_fig
#width/height: 1000/500