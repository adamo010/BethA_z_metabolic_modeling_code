#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 17:22:50 2021

@author: adamo010
"""

#based off of 07.09.21_Hale2018_OTU_table_filtering_pipeline_clean.py. NEW DATASET!
#based off of 10.22.21_Niccolai2020_OTU_table_filtering_pipeline.py. NEW DATASET AGAIN!

import scipy
import numpy as np
import csv
import subprocess
import os
import pandas as pd
import fnmatch
import matplotlib.pyplot as plt
import shutil
import glob
import re
import copy

#Sambhawa is removing them by only keeping taxa with relative abundances of >=0.001 in at least 20% of samples

############################Step 0: import and process the OTU table to a reasonable level############################
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal/")

preprocessed_otus = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/ASV_tables_Muehlbauer2021/11.30.21_DADA2_analysis_Muehlbauer_data_microbiome_ASVs_counts.tsv", sep="\s+") 
#this is so dumb, this file is space separated

#need to edit header names
preprocessed_otus.rename(columns=lambda x: x.rstrip("_F_trimmed_filtered.fastq.gz"), inplace=True)

#phew, okay, need to import taxonomy here also and merge it together.
taxonomy = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/ASV_tables_Muehlbauer2021/11.30.21_DADA2_analysis_Muehlbauer_data_microbiome_ASVs_taxonomy.tsv", sep="\s+")
#nice

#MERGE!
taxsplit_otus = pd.merge(right= taxonomy, left=preprocessed_otus, how= "left", left_index = True, right_index=True) #merging on indices here, not col names

#taxsplit_otus = copy.deepcopy(preprocessed_otus)   #no need for this anymore
taxsplit_otus.columns = taxsplit_otus.columns.str.replace(' ', '') #remove spaces from column headers
taxsplit_otus= taxsplit_otus.rename(columns={"domain": "kingdom"}) #rename columns (rename domain for ease of downstream analyses)

#Now, to create columns where I stitch these back together. Preferrably in function form.
#with the 03.09.20 addition of the blanked cells, I need to figure out some sort of conditional, where stuff is only added if 
#the cell isn't blank.
#basically, if taxsplit_otus (next) is blank, don't add anything. Otherwise, do the command

taxsplit_otus["phylum_level_sort"] = np.where(len(taxsplit_otus["phylum"]) > 1, (taxsplit_otus["kingdom"] + "; " + taxsplit_otus['phylum']), "")
taxsplit_otus["class_level_sort"] = np.where(len(taxsplit_otus["class"]) > 1, (taxsplit_otus["phylum_level_sort"] + "; " + taxsplit_otus['class']), "")
taxsplit_otus["order_level_sort"] = np.where(len(taxsplit_otus["order"]) > 1, (taxsplit_otus["class_level_sort"] + "; " + taxsplit_otus['order']), "")
taxsplit_otus["family_level_sort"] = np.where(len(taxsplit_otus["family"]) > 1, (taxsplit_otus["order_level_sort"] + "; " + taxsplit_otus['family']), "")
taxsplit_otus["genus_level_sort"] = np.where(len(taxsplit_otus["genus"]) > 1, (taxsplit_otus["family_level_sort"] + "; " + taxsplit_otus['genus']), "")
taxsplit_otus["species_level_sort"] = np.where(len(taxsplit_otus["species"]) > 1, (taxsplit_otus["genus_level_sort"] + ";" + taxsplit_otus['species']), "")
#note that we're keeping species_level_sort for organizing purposes only
#this dataset does NOT have a species column.

#Now I have an OTU table that I can use for sorting, filtering, etc. I will use this moving forward

############################Step 1: remove all obvious contaminants############################
#here is my contaminants list:
contaminants = pd.read_csv("06.16.20_Burns_contaminant_list.txt", sep="/t", header = None, engine='python') #new contaminants list
contaminants[0] = contaminants[0].str.strip()       #added 04.15.20
#specify the engine because an error message told me to 
#note that the contaminants list has to be a single column, in the same format as an OTU table spit out from QIIME (or, basically, the same format as preprocessed_otus)

#do this based on xx_level_sort- first, split 'contaminants' into different lists
kingdom_contaminants = []
phylum_contaminants = []
class_contaminants = []
order_contaminants = []
family_contaminants = []    #NEW from 04.27.20
genus_contaminants = []     #NEW from 04.27.20
for elem in contaminants[0]:
    if "g__" in elem:
        genus_contaminants.append(elem)#NEW from 04.27.20
    elif "f__" in elem:
        family_contaminants.append(elem)#NEW from 04.27.20
    elif "o__" in elem:
        order_contaminants.append(elem)
    elif "c__" in elem:
        class_contaminants.append(elem)
    elif "p__" in elem:
        phylum_contaminants.append(elem)
    elif "k__" in elem:
        kingdom_contaminants.append(elem)
    del elem

#NEW as of 07.09.21: letters and underscores from taxonomic ids. Not needed for Hale2018 taxonomy data. 
kingdom_contaminants = [s.replace("k__", "") for s in kingdom_contaminants]
phylum_contaminants = [s.replace("k__", "").replace("p__","").replace("[","").replace("]","") for s in phylum_contaminants]
class_contaminants = [s.replace("k__", "").replace("p__","").replace("c__","").replace("[","").replace("]","") for s in class_contaminants]
order_contaminants = [s.replace("k__", "").replace("p__","").replace("c__","").replace("o__","").replace("[","").replace("]","") for s in order_contaminants]
family_contaminants = [s.replace("k__", "").replace("p__","").replace("c__","").replace("o__","").replace("f__","").replace("[","").replace("]","") for s in family_contaminants]
genus_contaminants = [s.replace("k__", "").replace("p__","").replace("c__","").replace("o__","").replace("f__","").replace("g__","").replace("[","").replace("]","") for s in genus_contaminants]

#NEW as of 10.22.21: add an additional phylum contaminant, Bacteria; Cyanobacteria/Chloroplast;, that is specific to Niccolai2020 dataset but NOT in Muehlbauer
#phylum_contaminants.append("Bacteria; Cyanobacteria/Chloroplast") 

#now, run taxsplit_otus by each of these lists and remove all rows that match our contaminants list.
#at the same time, pull out contaminants to keep track of OTUs
taxsplit_otus_filtered = copy.deepcopy(taxsplit_otus) #create copy first
taxsplit_otus_contaminants= copy.deepcopy(taxsplit_otus)

#kingdom contaminants
taxsplit_otus_filtered= taxsplit_otus_filtered[~taxsplit_otus_filtered["kingdom"].isin(kingdom_contaminants)] #non-contaminants
taxsplit_otus_contaminants_k = taxsplit_otus_contaminants[taxsplit_otus_contaminants["kingdom"].isin(kingdom_contaminants)]

#phylum contaminants
taxsplit_otus_filtered= taxsplit_otus_filtered[~taxsplit_otus_filtered["phylum_level_sort"].isin(phylum_contaminants)] 
taxsplit_otus_contaminants_p = taxsplit_otus_contaminants[taxsplit_otus_contaminants["phylum_level_sort"].isin(phylum_contaminants)]

#class contaminants
taxsplit_otus_filtered= taxsplit_otus_filtered[~taxsplit_otus_filtered["class_level_sort"].isin(class_contaminants)] 
taxsplit_otus_contaminants_c = taxsplit_otus_contaminants[taxsplit_otus_contaminants["class_level_sort"].isin(class_contaminants)]

#order contaminants
taxsplit_otus_filtered= taxsplit_otus_filtered[~taxsplit_otus_filtered["order_level_sort"].isin(order_contaminants)] 
taxsplit_otus_contaminants_o = taxsplit_otus_contaminants[taxsplit_otus_contaminants["order_level_sort"].isin(order_contaminants)]

#family contaminants
taxsplit_otus_filtered= taxsplit_otus_filtered[~taxsplit_otus_filtered["family_level_sort"].isin(family_contaminants)] 
taxsplit_otus_contaminants_f = taxsplit_otus_contaminants[taxsplit_otus_contaminants["family_level_sort"].isin(family_contaminants)]

#genus contaminants
taxsplit_otus_filtered= taxsplit_otus_filtered[~taxsplit_otus_filtered["genus_level_sort"].isin(genus_contaminants)]
taxsplit_otus_contaminants_g = taxsplit_otus_contaminants[taxsplit_otus_contaminants["genus_level_sort"].isin(genus_contaminants)]

#track up all the contaminants by appending them into a single dataframe.
pd_list = [taxsplit_otus_contaminants_k, taxsplit_otus_contaminants_p, taxsplit_otus_contaminants_c, taxsplit_otus_contaminants_o,
           taxsplit_otus_contaminants_f,taxsplit_otus_contaminants_g]
taxsplit_otus_contaminants_all = pd.concat(pd_list)

#save both of these.
taxsplit_otus_filtered.to_csv("12.01.21_Muehlbauer2021_OTUs_contaminants_removed.csv")
taxsplit_otus_contaminants_all.to_csv("12.01.21_Muehlbauer2021_OTUs_contaminant_taxa.csv")

#NEW 07.13.21: adding a command to create a taxonomy table for lookups
taxsplit_otus["taxonomy_orig"] = taxsplit_otus['kingdom'].astype(str)+str("; ")+taxsplit_otus['phylum'].astype(str)+str("; ")\
    +taxsplit_otus['class'].astype(str)+str("; ")+taxsplit_otus['order'].astype(str)+str("; ")+taxsplit_otus['family'].astype(str)\
        +str("; ")+taxsplit_otus['genus'].astype(str)+str("; ")+taxsplit_otus['species'].astype(str)
#save this as a general reference in published_data
taxsplit_otus.to_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/ASV_tables_Muehlbauer2021/Muehlbauer2021_OTUs_with_matched_taxonomy.csv")

#cleanup:
del preprocessed_otus, kingdom_contaminants, phylum_contaminants, class_contaminants, order_contaminants, family_contaminants, genus_contaminants, taxsplit_otus, contaminants
del taxsplit_otus_contaminants_all, taxsplit_otus_contaminants_c, taxsplit_otus_contaminants_f, taxsplit_otus_contaminants_g, taxsplit_otus_contaminants_k,\
    taxsplit_otus_contaminants_p, taxsplit_otus_contaminants_o

############################Step 2: add in read counts############################
    #did this above, for Muehlbauer data, which used dada2 output, so we'll go with that. 

#step 0: import read counts table, remove contaminant taxa, merge with taxonomy, and move forward.

#read_counts = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Hale2018/sOTU_read_counts_by_sample_for_MICOM.csv") 
#read_counts= read_counts.drop(columns=["Unnamed: 0"])
#read_counts.set_index("OTU_ID_reads", inplace=True)
#read_counts=read_counts.T.reset_index().rename(columns={"index": "OTU_ID"}) #transpose, reset index, and rename first column
#contaminant_otus_list = taxsplit_otus_contaminants["OTU_ID"].to_list() #get list of OTUids identified as contaminants
#read_counts_filtered = copy.deepcopy(read_counts)
#read_counts_filtered = read_counts_filtered[~read_counts_filtered["OTU_ID"].isin(contaminant_otus_list)]

#########KEEP: for adding taxonomy data
#otus_plus_tax_filtered = pd.merge(left= read_counts_filtered, right=taxsplit_otus_filtered, how= "left", left_on ="OTU_ID", right_on="OTU_ID")
#otus_plus_tax_filtered.set_index("OTU_ID", inplace=True)
#add a "taxonomy_orig" column for collapsing whole taxonomy in subsequent data
taxsplit_otus_filtered["taxonomy_orig"] = taxsplit_otus_filtered['kingdom'].astype(str)+str("; ")+taxsplit_otus_filtered['phylum'].astype(str)\
    +str("; ")+taxsplit_otus_filtered['class'].astype(str)+str("; ")+taxsplit_otus_filtered['order'].astype(str)+str("; ")\
        +taxsplit_otus_filtered['family'].astype(str)+str("; ")+taxsplit_otus_filtered['genus'].astype(str)+str("; ")+\
            taxsplit_otus_filtered['species'].astype(str)

#cleanup
del taxsplit_otus_contaminants, pd_list

############################Step 3: abundance filtering############################
#step 0: get a list of sample names
otus_plus_tax_filtered = copy.deepcopy(taxsplit_otus_filtered)
otus_plus_tax_filtered= otus_plus_tax_filtered.reset_index()
sample_names = []
for col in otus_plus_tax_filtered.columns:       #use raw read counts table
    sample_names.append(str(col))
    del col
del sample_names[0:1] #delete first element, which is index
del sample_names[-14:] #delete non-sample name headers like genus (the cost of using sample names)
#for this dataset, should end up with a list of 12 sampleids (SRR####)'
#NOTE here- I'm also going to edit the column name of the dataframe. These data use ASVs, but... we're using OTU because I said so.
otus_plus_tax_filtered.rename(columns={"index":"OTU_ID"}, inplace=True)

#step 1: collapse by OTU_id
otu_table_collapsed = otus_plus_tax_filtered.groupby(['taxonomy_orig', 'OTU_ID'])[sample_names].sum()

#step 2: calculate relative abundances
otu_table_relabund = otu_table_collapsed.loc[:,].div(otu_table_collapsed.sum(axis=0))
#it is worth noting that this works BECAUSE in the groupby function in step 1, taxonomy_orig and #OTU_ID are both

#step 3: filter by abundance (at least 0.001 in at least 10% of samples)
num_above_cutoff = [(otu_table_relabund[sample_names] >= 0.001).sum(1)]
#convert to a list
for elem in num_above_cutoff:
    otu_table_num_above_cutoff = elem.values.tolist()       #this converts a series to a list
del elem
#add this list as a column
otu_table_relabund['Num_above_cutoff'] = otu_table_num_above_cutoff

#step 4: calculate the number of samples in the total dataframe- will need this value later
#use OTU_table_collapsed; for this dataframe, the sample names don't have a consistent character (like Sample_ in Burns data) so just
#count the number of columns in otu_table_collapsed, which only has samples as columns (the rest is set as index)
numsamples =0
for col in otu_table_collapsed:
    numsamples = numsamples+1
del col

#step 5: apply the abundance filter to the OTU table
otu_table_relabund_above_cutoff = otu_table_relabund.loc[otu_table_relabund['Num_above_cutoff'] >= (numsamples*0.1)]
#NEW as of 06.16.20: find out how many samples were abundance-filtered out. 
otu_table_relabund_below_cutoff = otu_table_relabund.loc[otu_table_relabund['Num_above_cutoff'] < (numsamples*0.1)]
#above cutoff= relabund3; belowcutoff= relabund4

#reset the indices
otu_table_relabund_above_cutoff.reset_index(inplace=True)
otu_table_relabund_below_cutoff.reset_index(inplace=True)
#NEW as of 06.16.20: find out how many samples were abundance-filtered out. 

#step 6: re-import the phylogeny stuff; BE CAREFUL about matching index with taxsplit_otus_filtered
otu_table_relabund_above_cutoff_plus_taxonomy = pd.merge(left= otu_table_relabund_above_cutoff, right=taxsplit_otus_filtered, how= "left", right_index=True, left_on="OTU_ID")
otu_table_relabund_below_cutoff_plus_taxonomy = pd.merge(left= otu_table_relabund_below_cutoff, right=taxsplit_otus_filtered, how= "left", right_index=True, left_on="OTU_ID")
#whew, it's a lot easier to do it this way than I did it in the 06.16.20 version, huh?

#step 7: save abundance data
otu_table_relabund_above_cutoff_plus_taxonomy.to_csv("12.01.21_Muehlbauer2021_asv_table_abundance_filtered_0.1.csv")
otu_table_relabund_below_cutoff.to_csv("12.01.21_Muehlbauer2021_asv_table_below_0.1_abundance_taxa.csv")

#step 8 NEW!!: filter by abundance round 2 (at least 0.001 in at least one sample)
#first, identify rows (OTUs) that do NOT hit at least 0.001 for every sample
otu_table_relabund2 = otu_table_collapsed.loc[:,].div(otu_table_collapsed.sum(axis=0))
keep_or_toss = []
for index, row in otu_table_relabund2.iterrows():
    rowlist = row.to_list() #convert each row to list
    if any(elem >= 0.001 for elem in rowlist):
        keep_or_toss.append("KEEP")
    else:
        keep_or_toss.append("TOSS") #scan through row lists and identify if any elemeng is >=0.001
    #print(rowlist)
otu_table_relabund2["keep_or_toss"] = keep_or_toss     #append this list to relabund table

#filter df by keep_or_toss list (i.e. abundance filter)
otu_table_relabund_filtered = otu_table_relabund2[otu_table_relabund2['keep_or_toss'].str.contains("KEEP")] 
otu_table_relabund_filtered_lowabundance = otu_table_relabund2[otu_table_relabund2['keep_or_toss'].str.contains("TOSS")] #NEW to track lowab asvs

#save 
otu_table_relabund_filtered.to_csv("12.01.21_Muehlbauer2021_asv_table_abundance_filtered_0.1_any_sample_cutoff.csv")
otu_table_relabund_filtered_lowabundance.to_csv("12.01.21_Muehlbauer2021_asv_table_abundance_below_0.1_any_sample_cutoff.csv")


############################################################################################################################################
#########################################Relative abundance filtering by taxonomic level########################################################
############################################################################################################################################

############################Step 1: create separate dataframes to filter/collapse at each taxonomic level############################
#get a list of the column names and their indices- will need these (probably) to make new dataframes
col_names_list = []
for col in otus_plus_tax_filtered.columns:
    col_names_list.append(str(col))
    del col
    
#create seperate dataframes 
kingdom_level_otus = otus_plus_tax_filtered.drop(["phylum", "class", "order", "family", "genus", "species", "phylum_level_sort", 
                                         "class_level_sort", "order_level_sort", "family_level_sort", "genus_level_sort", "species_level_sort"], axis = 'columns')    
phylum_level_otus = otus_plus_tax_filtered.drop(["kingdom", "class", "order", "family", "genus", "species",
                                        "class_level_sort", "order_level_sort", "family_level_sort", "genus_level_sort", "species_level_sort"], axis = 'columns')    
class_level_otus = otus_plus_tax_filtered.drop(["kingdom", "phylum", "order", "family", "genus", "species", "phylum_level_sort", 
                                        "order_level_sort", "family_level_sort", "genus_level_sort", "species_level_sort"], axis = 'columns')    
order_level_otus = otus_plus_tax_filtered.drop(["kingdom", "phylum", "class", "family", "genus", "species", "phylum_level_sort", 
                                       "class_level_sort", "family_level_sort", "genus_level_sort", "species_level_sort"], axis = 'columns')    
family_level_otus = otus_plus_tax_filtered.drop(["kingdom", "phylum", "class", "order", "genus", "species", "phylum_level_sort", 
                                        "class_level_sort", "order_level_sort", "genus_level_sort", "species_level_sort"], axis = 'columns')    
genus_level_otus = otus_plus_tax_filtered.drop(["kingdom", "phylum", "class", "order", "family", "species", "phylum_level_sort", 
                                       "class_level_sort", "order_level_sort", "family_level_sort", "species_level_sort"], axis = 'columns') 
species_level_otus = otus_plus_tax_filtered.drop(["kingdom", "phylum", "class", "order", "family", "genus", "phylum_level_sort", 
                                       "class_level_sort", "order_level_sort", "family_level_sort", "genus_level_sort"], axis = 'columns') 
   

############################Step 2: collapse at each taxonomic level############################

#note that kingdom level sorting only involves removing archaea- that was done in the contaminant removal step above. So kingdom code is short
kingdom_level_otus=kingdom_level_otus.sort_values(by=["kingdom"])

#step 0: pull all columns that need to be summed (i.e. all OTU columns)- will use a lot
otus_plus_tax_filtered = otus_plus_tax_filtered.reset_index()
#used to have a sample_names thing here, but I've generated it above, so I won't worry about it. 

def tax_level_filtering(tax_level, tax_level_otus, tax_level_sort):
    tax_level_collapsed = tax_level_otus.groupby([tax_level_sort])[sample_names].sum()
    tax_level_relabund = tax_level_collapsed.loc[:,].div(tax_level_collapsed.sum(axis = 0))
    num_above_cutoff = [(tax_level_relabund[sample_names] >= 0.001).sum(1)]
    for elem in num_above_cutoff:
        tax_level_num_above_cutoff = elem.values.tolist()
    tax_level_relabund["Num_above_cutoff"] = tax_level_num_above_cutoff
    numsamples = 0
    for col in tax_level_relabund:
        if 'Sample_' in str(col):
            numsamples = numsamples+1
    tax_level_relabund2 = tax_level_relabund.loc[tax_level_relabund['Num_above_cutoff'] >= (numsamples*0.1)]
    tax_level_relabund2.to_csv("_relabund_0.1_Muehlbauer2021.csv")  
    for file in os.listdir():                       #rename files by their appropriate diet; couldn't figure out a way to do this in the above code
        src=file
        if fnmatch.fnmatch(file, "_relabund_0.1_Muehlbauer2021.csv"):
            dst = str(tax_level)+str('_level')+file
            os.rename(src,dst)
    return

tax_level_filtering("phylum", phylum_level_otus, 'phylum_level_sort')
tax_level_filtering("class", class_level_otus, 'class_level_sort')
tax_level_filtering("order", order_level_otus, 'order_level_sort')
tax_level_filtering("family", family_level_otus, 'family_level_sort')
tax_level_filtering("genus", genus_level_otus, 'genus_level_sort')

#well, that's that. Much easier with this dataset, for some reason. 
