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

#starting point is 02.22.22_Burns2022_CORDA_model_analyzing_EUstd_uptakes

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_analyses/") 

##########Step 2: import data and adjust factors
upsecs_all <- read.csv(file= "03.21.22_combining_host_and_microbiomes_EUstd_host_uptakes_and_secretions.csv", header=TRUE)
upsecs_all <- subset(upsecs_all, select=-c(X))
upsecs_all$Patient_Blind_ID <- as.factor(upsecs_all$Patient_Blind_ID)
upsecs_all$host_model_id <- as.factor(upsecs_all$host_model_id)
upsecs_all$microbiome_model_id <- as.factor(upsecs_all$microbiome_model_id)

#I suspect we'll also want to create unique ids for all of our samples
#first, need to resort columns so the columns to be merged are next to each other
upsecs_all <- upsecs_all[, c(1,2,3,4,6,5,7,8)]
upsecs_all <- upsecs_all %>% unite("paired_id", host_model_id:microbiome_model_id, remove = FALSE)

#COUNT 1: how many unique fluxes are we including here? 
upsecs_unique <- unique(upsecs_all$metabolite_name) #123 uptake rxs
rm(upsecs_unique) #cleanup

#########Step 3: filter out zeroes
analyzed_upsecs_wide <- select(upsecs_all, -c(paired_id, host_model_id, microbiome_model_id, host_model_description, microbiome_model_description)) 
#first, convert to wide form, dropping all columns that will impact spreading to wide form
analyzed_upsecs_wide <- analyzed_upsecs_wide %>% 
  spread(combo_type, metab_amount)
#then, drop all rows where all four categories of values are nan.
analyzed_upsecs_wide_noNANs <- analyzed_upsecs_wide[!(is.na(analyzed_upsecs_wide$tumor_tumor) & 
                                                          is.na(analyzed_upsecs_wide$tumor_normal) &
                                                          is.na(analyzed_upsecs_wide$normal_tumor) &
                                                          is.na(analyzed_upsecs_wide$normal_normal)),]
#perhaps unsurprisingly, there are none. 
upsecs_unique_zerofree <- unique(analyzed_upsecs_wide_noNANs$metabolite_name) #count the number of metabs: 123
#for this dataset, there are no nans and all data are kept
#replace remaining NA values with zeroes
analyzed_upsecs_wide_noNANs$tumor_tumor[is.na(analyzed_upsecs_wide_noNANs$tumor_tumor)] = 0
analyzed_upsecs_wide_noNANs$normal_normal[is.na(analyzed_upsecs_wide_noNANs$normal_normal)] = 0
analyzed_upsecs_wide_noNANs$tumor_normal[is.na(analyzed_upsecs_wide_noNANs$tumor_normal)] = 0
analyzed_upsecs_wide_noNANs$normal_tumor[is.na(analyzed_upsecs_wide_noNANs$normal_tumor)] = 0

#########Step 4: filter by uptakes which are active in at least half the total number of samples 
#first, create a filtering dataframe that contains counts of each Reaction incidence in the dataframe
#need to have n=2 if values are nonzero in both tumor and normal columns; n=1 if values are nonzero in only one column
analyzed_upsecs_filtering_df <- analyzed_upsecs_wide_noNANs %>%  #create a new dataframe by filtering the old one
  mutate(count_of_nonzeros = rowSums(.[3:6]!=0)) %>% #create a new column, count_of_nonzeroes, that counts # nonzero values in columns 3 and 4 (tumor and normal) for each row
  group_by(metabolite_name) %>% #group by the column of interest
  summarize(n = sum(count_of_nonzeros)) %>% #sums the count_of_nonzeros values within each fluxpath_subsystem and saves the value as n
  mutate(n) #add n as a new column in the dataframe
#create a new dataframe that merges the old and new dataframes
analyzed_upsecs_filter_added <- dplyr::inner_join(analyzed_upsecs_wide_noNANs, analyzed_upsecs_filtering_df, by= "metabolite_name")
#OKAY one stupid thing here is that we have 44*4 = 176 possible values. So, I'm not exactly sure what the best filtering is here.
#n=176 means the metabolite is nonzero in every PBID in every host-microbiome combo.
#NEW step for combo data specifically; add an additional column that counts the number of combinations in each row
analyzed_upsecs_filter_added <- analyzed_upsecs_filter_added %>% mutate(count_of_nonzeros =rowSums(.[3:6]!=0)) 
#the goal is to have count_of_nonzeros = 4, which means the metabolite is nonzero in all four microbiome/host status

#okay, let's start filtering. 
#step1: remove all metabolites not persent in all four host/microbiome status samples
analyzed_upsecs_filtered1 <- filter(analyzed_upsecs_filter_added, count_of_nonzeros == 4)
#count number of unique metabolites left
upsecs_filtered1_unique <- unique(analyzed_upsecs_filtered1$metabolite_name) #this gets us to 44 uptake/secretions
#do a second filtering step: My guess is that we won't lose that many uptakes/secretions in this step
#we're goint to filter by 88 samples, which is half. 
analyzed_upsecs_filtered2 <- filter(analyzed_upsecs_filtered1, n >= 88)
#count number of unique metabolites left
upsecs_filtered2_unique <- unique(analyzed_upsecs_filtered2$metabolite_name) #that gets us to 24 metabolites. Probably fine. 

#clean up
rm(analyzed_upsecs_filtered1, upsecs_filtered1_unique, upsecs_filtered2_unique)

##########Step 5: split into multiple dataframes by metabolite and do stats on each metabolite
#first, edit the dataframe
#drop a couple of columns
analyzed_upsecs_paired_only <- select(analyzed_upsecs_filtered2, -c(count_of_nonzeros, n))
#rename a couple of columns to make life easier later. Yeah, yeah, I know, this is the reverse of what I did before.
#convert to long form
analyzed_upsecs_paired_long <- reshape2::melt(data= analyzed_upsecs_paired_only,
                                               id.vars= c("Patient_Blind_ID", "metabolite_name"),
                                               variable.name = "combo_type",
                                               value.name = "metabolite_amount")
#split up by different subsystems and see where we get. '
#NOTE here that we're running a Kruskal-wallis test; once we have some (hopefully) significant results, we'll do the pairwise stuff.
split_analyzed_upsecs_by_metab <- split(analyzed_upsecs_paired_long, with(analyzed_upsecs_paired_long, interaction(metabolite_name)), drop = TRUE)
by_metab_stats <- lapply(split_analyzed_upsecs_by_metab , function(df){kruskal.test(metabolite_amount ~ combo_type, data=df)})

#TOGGLE stats here

##########Step 6: statistics
#make lists of p values
by_metab_pvals <- c() #initialize list
for (elem in by_metab_stats){
  new_value = elem$p.value
  by_metab_pvals <- c(by_metab_pvals, new_value)}
rm(new_value)
#make a list of q-values
by_metab_qvals <- p.adjust(by_metab_pvals, method = "fdr")
#make a list of genus_ids
metab_ids <- c()
for(elem in split_analyzed_upsecs_by_metab){
  new_value = elem$metabolite_name[1]
  metab_ids <- c(metab_ids, new_value)}
#merge all lists together. 
metabolite_statistics <- data.frame(metab_ids, by_metab_pvals, by_metab_qvals)

#########Step 7: merge descriptions back in
metabolites_metadata_slimmed <- select(upsecs_all, c(metabolite_name, combo_type))
metabolites_metadata_slimmed <- metabolites_metadata_slimmed[!duplicated(metabolites_metadata_slimmed$metabolite_name), ]
metabolite_statistics_plus_metadata <-left_join(metabolite_statistics, metabolites_metadata_slimmed,
                                              by = c("metab_ids" = "metabolite_name"))
##########Step 8: save results as a csv file
metabolite_statistics_plus_metadata <- metabolite_statistics_plus_metadata[order(metabolite_statistics_plus_metadata$by_metab_qvals),] #re-order by q-value
metabolite_statistics_plus_metadata <- metabolite_statistics_plus_metadata %>% rename("metabolite_name"="metab_ids")
write.csv(as.data.frame(metabolite_statistics_plus_metadata), file="03.23.22_combining_host_and_microbiomes_EUstd_host_uptake_and_secretion_stats.csv")

##########Step 9: graph
#can just use analyzed_upsecs_paired_long, since we don't care about pairwise differences here
y_comparisons <- list( c("tumor_tumor", "tumor_normal"), c("tumor_tumor", "normal_tumor"), c("tumor_tumor", "normal_normal"),
                       c("tumor_normal", "normal_tumor"), c("tumor_normal", "normal_normal"), c("normal_tumor", "normal_normal"))
#graph: atp
ggplot(subset(analyzed_upsecs_paired_long, metabolite_name %in% c("EX_atp[e]")), aes(x= combo_type, y= metabolite_amount, color=combo_type)) +
  geom_violin(aes(x= combo_type, y= metabolite_amount, fill=combo_type)) + 
  geom_dotplot(binaxis='y', stackdir='center', stackratio=0.75, dotsize=0.6, aes(x= combo_type, y= metabolite_amount)) +
  labs(title = expression("Colonocyte ATP uptakes, Burns2015 data"), 
       subtitle= expression("p-values from uncorrected pairwise Mann-Whitney U tests")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Metabolite amounts,\nmmol/(gDW * h)", limits=c(-1000,800)) +
  scale_x_discrete(labels = c('Normal\nmicrobiome,\nnormal host','Normal\nmicrobiome,\ntumor host',
                              'Tumor\nmicrobiome,\nnormal host','Tumor\nmicrobiome,\ntumor host'),
                   name= "Patient microbiome and host tissue descriptions") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "none") +
  guides(fill='none') +
  stat_compare_means(comparisons= y_comparisons, method = "wilcox.test", size= 3)

#graph: utp
ggplot(subset(analyzed_upsecs_paired_long, metabolite_name %in% c("EX_utp[e]")), aes(x= combo_type, y= metabolite_amount, color=combo_type)) +
  geom_violin(aes(x= combo_type, y= metabolite_amount, fill=combo_type)) + 
  geom_dotplot(binaxis='y', stackdir='center', stackratio=0.75, dotsize=0.6, aes(x= combo_type, y= metabolite_amount)) +
  labs(title = expression("Colonocyte UTP uptakes, Burns2015 data"), 
       subtitle= expression("p-values from uncorrected pairwise Mann-Whitney U tests")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Metabolite amounts,\nmmol/(gDW * h)", limits=c(-1,4)) +
  scale_x_discrete(labels = c('Normal\nmicrobiome,\nnormal host','Normal\nmicrobiome,\ntumor host',
                              'Tumor\nmicrobiome,\nnormal host','Tumor\nmicrobiome,\ntumor host'),
                   name= "Patient microbiome and host tissue descriptions") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "none") +
  guides(fill='none') +
  stat_compare_means(comparisons= y_comparisons, method = "wilcox.test", size= 3, label.y = c(1.8,2.2,2.6,3.0,3.4,3.8))

#SUCKS
