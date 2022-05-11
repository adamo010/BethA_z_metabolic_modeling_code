library(dplyr)
library(tidyr)
library(ggplot2) # needs to be version ≥ 2.1.0
library(scales)
library(network)
library(sna)
library(GGally)
library(geomnet)
library(ggnetwork)
library(igraph)
library(dplyr)
library(tidyverse)
library(viridis)
library(ggpubr)

#borrowing graphing code from 04.20.22_ggnet2_network_plot_tutorial.R
#data for this script generated in 04.12.22_correlating_Burns_GEs_and_all_microbiome_metrics.R

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/") 

###########Import data
#recall that these are all correlations with UNCORRECTED p-values >0.01
gr_ge_corrs <- read.csv("04.12.22_Burns2015_metabolic_gene_expression_by_gr_sigto0.01_correlations.csv")
gr_ge_corrs <- subset(gr_ge_corrs, select = -c(X))
gr_ge_corrs <- gr_ge_corrs %>%                               
  mutate(taxon_name = stringr::str_split(taxon_name, "_") %>% map_chr(., 2)) #clip off Family names
flux_ge_corrs <- read.csv("04.12.22_Burns2015_metabolic_gene_expression_by_fluxpath_sigto0.01_correlations.csv")
flux_ge_corrs <- subset(flux_ge_corrs, select = -c(X))
#bring in some other fluxpath names
fluxpath_key <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Common_growth_medium_runs/May_June2021_MICOM_flux_and_pathway_analyses/Microbiome_flux_key_full_June2021.csv", sep=',')#there are a lot of extra columns here: ony really want the first 2
fluxpath_key <- subset(fluxpath_key, select= c("abbreviation", "description"))
flux_ge_corrs <- left_join(flux_ge_corrs, fluxpath_key, by= c("fluxpath_name" = "abbreviation"))
flux_ge_corrs <- flux_ge_corrs %>%                               
  mutate(metab = stringr::str_split(description, " of ") %>% map_chr(., 2)) #clean up metab names
flux_ge_corrs = select(flux_ge_corrs, metab, gene_id, p_value_adj, p_value_raw, R_value)

###########Create some network features
#convert node columns to character strings
gr_ge_corrs <- gr_ge_corrs %>%
  mutate_at(c("taxon_name", "gene_id"), as.character) #NOTE here we're using the microbial feature as the "from_id", and the host gene as the "to_id"

flux_ge_corrs <- flux_ge_corrs %>%
  mutate_at(c("metab", "gene_id"), as.character) #NOTE here we're using the microbial feature as the "from_id", and the host gene as the "to_id"

#first, we want to add a vertex attribute: color by host gene/microbe status
ge_gr_host <- unique(gr_ge_corrs$gene_id)
ge_gr_microbe <- unique(gr_ge_corrs$taxon_name)

ge_flux_host <- unique(flux_ge_corrs$gene_id)
ge_flux_microbe <- unique(flux_ge_corrs$metab)

#create a network
gr_ge.net <- network::network(gr_ge_corrs[, 1:2], directed = FALSE) #create an undirected network from the first two columns of dataframe

flux_ge.net <- network::network(flux_ge_corrs[, 1:2], directed = FALSE) #create an undirected network from the first two columns of dataframe

#add categorical vertex attribute
##creates a node_type column set as True if "network.vertex.names are found in ge_gr_host vector (a list of gene names)
gr_ge.net %v% "node_type" <-  (network.vertex.names(gr_ge.net) %in% ge_gr_host) 
gr_ge.net %v% "node_type" <-  1 + as.integer(gr_ge.net %v% "node_type") #converts True to 2, and False to 1
rownames(gr_ge_corrs$node_type) <- gr_ge_corrs$node_type$name

flux_ge.net %v% "node_type" <-  (network.vertex.names(flux_ge.net) %in% ge_flux_host) 
flux_ge.net %v% "node_type" <-  1 + as.integer(flux_ge.net %v% "node_type") #converts True to 2, and False to 1
rownames(flux_ge_corrs$node_type) <- flux_ge_corrs$node_type$name

#double check; should have new column called node_type where genes are marked as 2, and bacteria are marked as 1
head(ggnetwork(gr_ge.net))
tail(ggnetwork(gr_ge.net))

head(ggnetwork(flux_ge.net))
tail(ggnetwork(flux_ge.net))

## create edge attribute (R value)
network::set.edge.attribute(gr_ge.net, "R_value", gr_ge_corrs$R_value)

network::set.edge.attribute(flux_ge.net, "R_value", flux_ge_corrs$R_value)

#all right. Now, what I want to do is import the actual correlation data and add an edge attribute that indicates the correlation direction.
#lucky for me, the R-value gives that information. add a new column to reflect that. 
gr_ge_corrs <- gr_ge_corrs %>%
  mutate(corr_dir = case_when(
    R_value <= 0 ~ "negative",
    R_value > 0 ~ "positive"
  ))

flux_ge_corrs <- flux_ge_corrs %>%
  mutate(corr_dir = case_when(
    R_value <= 0 ~ "negative",
    R_value > 0 ~ "positive"
  ))

#now, add new edge attribute
network::set.edge.attribute(gr_ge.net, "corr_dir", gr_ge_corrs$corr_dir)
network::set.edge.attribute(flux_ge.net, "corr_dir", flux_ge_corrs$corr_dir)

#let's graph again
#optimal size is 700x700
plot1 <- ggplot(data= ggnetwork(gr_ge.net, layout.alg = "eigen$symstrong"), 
                aes(x, y, xend=xend, yend=yend)) +
  geom_edges(aes(size = R_value, color= corr_dir)) + #sets edge widths to n
  geom_nodes(aes(colour=factor(node_type)), size = 5) + 
  geom_nodelabel_repel(aes(label= vertex.names, colour=factor(node_type), size=2.5, fontface= "italic"), 
                       force = 5, label.padding=unit(0.15, "cm"), point.padding = unit(0.15, "cm")) +
  theme_blank() +
  theme(legend.position = "none", plot.margin = margin(1,1,1,1, "cm")) +
  scale_fill_manual(values=c("#332288","#117733","#AA4499","#E69F00")) + #bacteria, host, negative, positive
  scale_color_manual(values=c("#332288","#117733","#AA4499","#E69F00")) 

plot2 <- ggplot(data= ggnetwork(flux_ge.net, layout.alg = "eigen$symstrong"), #spring, rmds, fruchtermanreingold
                aes(x, y, xend=xend, yend=yend)) +
  geom_edges(aes(size = R_value, color= corr_dir)) + #sets edge widths to n
  geom_nodes(aes(colour=factor(node_type)), size = 5) + 
  geom_nodelabel_repel(aes(label=vertex.names, colour=factor(node_type), size=2), fontface= "italic",
                       force = 3, label.padding=unit(0.15, "cm"), point.padding = unit(0.5, "cm")) +
  theme_blank() +
  theme(legend.position = "none", plot.margin = margin(1,1,1,1, "cm")) +
  scale_fill_manual(values=c("#332288","#117733","#AA4499","#E69F00")) + #bacteria, host, negative, positive
  scale_color_manual(values=c("#332288","#117733","#AA4499","#E69F00")) 

#graph into a figure
ggarrange(plot1, plot2, 
          labels = c("A.", "B."),
          ncol = 2, nrow = 1)
#optimal size 1350, 700



####################Other color options
  scale_fill_manual(values=c("#61D04F","#2297E6","#DF536B","#F5C710")) + #bacteria, host, negative, positive
  scale_color_manual(values=c("#61D04F","#2297E6","#DF536B","#F5C710")) 

  scale_fill_manual(values=c("#10a53dFF","#3a5e8cFF","#ffcf20ff","#541352FF")) + #bacteria, host, negative, positive
    scale_color_manual(values=c("#10a53dFF","#3a5e8cFF","#ffcf20ff","#541352FF")) 

  scale_fill_manual(values=c("#e37e00","#1a476f","#c10534","#55752f")) + #bacteria, host, negative, positive
    scale_color_manual(values=c("#e37e00","#1a476f","#c10534","#55752f")) 
  
  scale_fill_manual(values=c("#aa4499","#332288","#d55e00","#117733")) + #bacteria, host, negative, positive
    scale_color_manual(values=c("#aa4499","#332288","#d55e00","#117733")) 
 
  scale_fill_manual(values=c("#785ef0","#dc2671","#fe6100","#648fff")) + #bacteria, host, negative, positive
    scale_color_manual(values=c("#785ef0","#dc2671","#fe6100","#648fff")) 

  scale_fill_manual(values=c("#332288","#117733","#AA4499","#E69F00")) + #bacteria, host, negative, positive
    scale_color_manual(values=c("#332288","#117733","#AA4499","#E69F00")) 
  
  
  #scale_fill_viridis(discrete= TRUE) +
  #scale_color_viridis(discrete= TRUE) 
  
  