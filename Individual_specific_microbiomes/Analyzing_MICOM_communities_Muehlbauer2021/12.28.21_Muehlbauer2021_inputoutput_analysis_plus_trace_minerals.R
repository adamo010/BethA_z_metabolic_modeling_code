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

#the goal here is to calculate p-values for the significant differences between input/output metabolites to use for correlating with gene expression
#differences. Also, we want to abundance filter input/output metabolites that do not appear in at least half the number of samples 
#(before running stats on identifying metabolites with significantly different input/outputs between tumor and normal samples).

#loosely based on 05.14.21_microbial_GR_analysis_Burns_data_filtered_top_20_taxa.R

#further based on Date_TBA_microbial_inputoutput_analysis_Burns_data_filtered_top_20_metabolites.R, with updates for Hale data
#and 09.07.21_Hale2018_inputoutput_analysis.R
#and 09.23.21_Hale2018_inputoutput_analysis_plus_trace_minerals.R

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Niccolai2020/") 

##########Step 1: import data and adjust factors
#It is important to note here that I have data (and have done analyses) for both family and genus level. However, I went with genus level
#moving forward because... most of the family level associations are basically due to a single genus anyway, I think is what I decided. 
analyzed_metabolites <- read_csv("11.09.21_Niccolai2020_inputoutput_with_metadata_longform.csv")
analyzed_metabolites$Patient_Blind_ID <- as.factor(analyzed_metabolites$Patient_Blind_ID)
#colnames(analyzed_metabolites) #print column names as needed

#how many unique metabolites are we including here? As of 09.22.21, 130 metabolites
metabolites_unique <- unique(analyzed_metabolites$metabolite_name)

#########Step 2 NEW as of 09.20.21: collapsing by trace minerals
#there are several metabolites (Calcium, Chloride ion, Cobalt, Cu2+, Heme, Magnesium, Mn2+, Potassium, Siroheme, and Zinc (II) ion)
#which have identical fluxes in all communities. The goal is to collapse these into a single group, called "Trace minerals"
#first, create a vector of these metabolites
trace_minerals <- c("Calcium", "Chloride ion", "Cobalt", "Cu2+", "Heme", "Magnesium", "Mn2+", "Potassium", "Siroheme", "Zinc (II) ion")
#Figure out some way to count the number of trace minerals with nonzero values within each sample. 
#first, create a dataframe that ONLY contains trace minerals
analyzed_metabolites_trace_minerals_only <- analyzed_metabolites[analyzed_metabolites$metabolite_name %in% trace_minerals, ]

#Use dplyr to round to 6 significant digits. That seems good enough
analyzed_metabolites_trace_minerals_only <- analyzed_metabolites_trace_minerals_only %>% 
  mutate(metabolite_amount=signif(metabolite_amount, 6)) 

#rerun the filtering_DF
trace_minerals_filtering_df <- analyzed_metabolites_trace_minerals_only %>%  #create a new dataframe by filtering the old one
  group_by(Sample_ID) %>% #group by the column of interest
  summarize(n = n_distinct(metabolite_amount)) #count the number of occurrences of each unique value in metabolite_amount and store as variable 'n'

#okay, now we have all 1s, with one 3 that is... the differences are too small to be relevant.So... let's just ignore it.

#time to collapse! for each sample_ID_ish, only keep ONE of the rows where metabolite_names values are in the trace_minerals vector
#and rename metabolite_name to trace_mineral.
analyzed_metabolites_collapsed <-data.frame(analyzed_metabolites) #copy the original dataframe
#first, rename "Potassium" to "trace elements"
analyzed_metabolites_collapsed$metabolite_name[analyzed_metabolites_collapsed$metabolite_name== "Potassium"] <- "Trace_elements"
#then, remove all rows where trace_minerals vector values are in the metabolite_name column
analyzed_metabolites_collapsed<- analyzed_metabolites_collapsed[!(analyzed_metabolites_collapsed$metabolite_name %in% trace_minerals),]

#NICE! Let's move forward. 
#how many unique metabolites are we including here? As of 09.22.21, 130 metabolites
metabolites_collapsed_unique <- unique(analyzed_metabolites_collapsed$metabolite_name)

#########Step 2: filter to only include paired taxa (i.e. taxa which appear in both tumor and normal samples)- used to be step 3
analyzed_metabolites_wide <- select(analyzed_metabolites_collapsed, -c(...1, Sample_ID)) #drop unnecessary columns
#first, convert to wide form
analyzed_metabolites_wide <- analyzed_metabolites_wide %>% 
  spread(Description, metabolite_amount)
#then, drop all rows where both tumor and normal values are nan.
analyzed_metabolites_wide_noNANs <- analyzed_metabolites_wide[!(is.na(analyzed_metabolites_wide$tumor) & is.na(analyzed_metabolites_wide$normal)),]
metabolites_unique_zerofree <- unique(analyzed_metabolites_wide_noNANs$metabolite_name) #count the number of metabolites

#replace remaining NA values with zeroes
analyzed_metabolites_wide_noNANs$tumor[is.na(analyzed_metabolites_wide_noNANs$tumor)] = 0
analyzed_metabolites_wide_noNANs$normal[is.na(analyzed_metabolites_wide_noNANs$normal)] = 0

#create a difference column
analyzed_metabolites_wide_noNANs$difference <- (analyzed_metabolites_wide_noNANs$tumor - analyzed_metabolites_wide_noNANs$normal)

#remove any rows where differences are 0
analyzed_metabolites_wide_noNANs_zerofree <- analyzed_metabolites_wide_noNANs[(analyzed_metabolites_wide_noNANs$difference !=0),]
#NOTE: for this dataset, nothing is removed.

#########Step 3: filter by metabolites which appear in at least half the total number of samples 
#we'll set a cutoff for at least 10-6.
#Filter to only include metabolites which appear in more than half of the samples, which translates to... 
#number of times a value appears greater than 10-6 in metabolite_name 
#first, create a filtering dataframe that contains counts of each metabolite_name incidence in the dataframe
analyzed_metabolites_filtering_df <- analyzed_metabolites_wide_noNANs_zerofree %>%  #create a new dataframe by filtering the old one
  mutate(count_of_nonzeros = rowSums(.[3:4]!=0)) %>% #create a new column, count_of_nonzeroes, that counts # nonzero values in tumor and normal flux columns
  group_by(metabolite_name) %>% #group by the column of interest
  summarize(n = sum(count_of_nonzeros)) %>% #sums the count_of_nonzeros values within each metabolite_name and saves the value as n
  mutate(n) #add n as a new column in the dataframe

#create a new dataframe that merges the old and new dataframes
analyzed_metabolites_filter_added <- dplyr::inner_join(analyzed_metabolites_wide_noNANs_zerofree, analyzed_metabolites_filtering_df, by= "metabolite_name")

#now, remove all rows where n (the name of the count column) is less than 16 (20% of the total number of samples)
analyzed_metabolites_filtered <- filter(analyzed_metabolites_filter_added, n >= 16)
#count number of unique metabolites left
metabolites_filtered_unique <- unique(analyzed_metabolites_filtered$metabolite_name)

##########Step 4: split into multiple dataframes by metabolite and do stats on each metabolite
#first, edit the dataframe
#drop a couple of columns
analyzed_metabolites_paired_only <- select(analyzed_metabolites_filtered, -c(difference, n))
#rename a couple of columns to make life easier later. Yeah, yeah, I know, this is the reverse of what I did before.
#analyzed_metabolites_paired_only <- analyzed_metabolites_paired_only %>% rename(normal= mean_genus_GR_normal, tumor= mean_genus_GR_tumor) #rename these columns
#convert to long form
analyzed_metabolites_paired_long <- reshape2::melt(data= analyzed_metabolites_paired_only,
                                              id.vars= c("Patient_Blind_ID", "metabolite_name"),
                                              variable.name = "Description",
                                              value.name = "metabolite_amount")
#split up by different metabolites and see where we get. 
split_analyzed_metabolites_by_metab <- split(analyzed_metabolites_paired_long, with(analyzed_metabolites_paired_long, interaction(metabolite_name)), drop = TRUE)
by_metab_stats <- lapply(split_analyzed_metabolites_by_metab, function(df){wilcox.test(metabolite_amount~Description, data=df, exact= FALSE, paired= TRUE)})

##########Step 5: statistics
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
for(elem in split_analyzed_metabolites_by_metab){
  new_value = elem$metabolite_name[1]
  metab_ids <- c(metab_ids, new_value)}
#merge all lists together. 
metabolite_statistics <- data.frame(metab_ids, by_metab_pvals, by_metab_qvals)

###########Step 6: clean up and remove extra variables
rm(elem, split_analyzed_metabolites_by_metab, 
   by_metab_pvals, by_metab_qvals, metab_ids, new_value)

###########Step 7: pull out significant p-values only
metabolite_statistics <- metabolite_statistics[order(metabolite_statistics$by_metab_pvals),] #re-order by p-value
#keep all 22. Whatever. 

##########Step 8: save results as a csv file
write.csv(as.data.frame(metabolite_statistics), file="11.09.21_filtered_microbial_input_output_all_results_Niccolai2020_data.csv")
#new: filter by sig p-values and only keep significantly different p-values
#metabolite_statistics_sig <- metabolite_statistics[metabolite_statistics$by_metab_pvals <= 0.05, ]
#write.csv(as.data.frame(metabolite_statistics_sig), file="11.09.21_filtered_microbial_input_output_SIG_results_Hale2018_data.csv")

##########Step 9: prepare data for graphing 
#LAAAAAAAAAAAAME now I have to calculate differences.
#use analyzed_metabolites_paired_only. First, drop boring cols
#analyzed_metabolites_paired_wide = subset(analyzed_metabolites_paired_only, select = -c(n,SampleID) )
#Split up to wide.
analyzed_metabolites_paired_wide <- spread(analyzed_metabolites_paired_long, Description, metabolite_amount)
#add a new column called diff. tumor-normal
analyzed_metabolites_paired_wide$diff <- (analyzed_metabolites_paired_wide$tumor - analyzed_metabolites_paired_wide$normal)

#NOW we can add a new column with useful binning of differences for graphing purposes
analyzed_metabolites_paired_wide <- analyzed_metabolites_paired_wide %>%
  mutate(Difference_direction = case_when(
    diff > 0 ~ "tumor > normal",
    diff < 0 ~ "tumor < normal",
    diff == 0 ~ "Zero change",
  ))

write.csv(as.data.frame(analyzed_metabolites_paired_wide), file="11.09.21_filtered_microbial_input_output_all_data_Niccolai2020_data.csv")

#step 10: NOW GRAPH??!!!
#count number of samples with Fuso:
freqs <- table(analyzed_metabolites_paired_wide$metabolite_name)
freqs

#Trace minerals
ggpaired(subset(analyzed_metabolites_paired_wide, metabolite_name %in% c("Trace_elements")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Trace elements exchange, Niccolai2020 data"), subtitle = expression("paired Wilcoxan, p=0.0158, q= 0.127")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(0,0.0008)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

