#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 15:27:57 2020

@author: adamo010
"""
#NEW as of 06.16.20: find out how many samples were abundance-filtered out. ALSO, new contaminants table.

#edits from 04.15.20 version: found some contaminant OTUs that slipped through the cracks (probably due to some weird spacing issues)
#some of these OTUs are below the order level; some are family, some are genus; probably need to add code to fix these

#edits from 04.13.20 version: onlu ~100 OTUs now have models assigned. Trying to figure out why. 

#edits from 03.09.20 version: ran THE WHOLE MODEL PICKING PIPELINE TWICE and found out that exactly the issue I thought
#would arise... did arise. Basically, a bunch of stuff came to the surface because the contaminants list I built previously
#was built on an already abundance-filtered table. I also moved files around so if there's an error, that's probably why.
#the desired filepath is now /Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal 

#edits from 03.04.20 version: add a 'removal of c__, f__ etc.' step. Needed later so should be done here.

###Edits from 02.24.20 version to match Sambhawa's approach: remove archea, then filter by abundance, then remove problem taxa
###Part of me doesn't really want to do this- I don't know if filtering by abundance before removing the problem taxa is the right
###way to do this. That's how you get weird taxa reappearing.

#compare my 02.24.20 list to Sambhawa's 02.26.20 list (from lab meeting)
#conclusions in 03.04.20_difficult_contaminant_taxa_Burns2015.xlsx
#Comparison in Burns_contams_removed_from_02.21_contam_list.xlsx
#came up with a new list. contaminants should now refer to 03.04.20_Burns_CRC_contaminants_list.txt

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
preprocessed_otus = pd.read_csv("Burns_absolute_OTU_table.csv") #note that this is in /Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal

taxsplit_otus = copy.deepcopy(preprocessed_otus)   #copy for safety
taxsplit_otus['taxonomy_orig'] = taxsplit_otus['taxonomy']     #copy the taxonomy column to save it for later
taxsplit_otus["taxonomy"] = taxsplit_otus["taxonomy"].str.split(";", expand=False)    #split the values in the taxonomy column into a list; delimiter is a semicolon

#now, to create new columns based on different taxonomic levels
def taxonomy_level_column_creation(otu_table, elem_num, tax_level_name):
    tax_level_list = []
    for elem in otu_table["taxonomy"]:
        tax_level_list.append(elem[elem_num])
    otu_table[str(tax_level_name)] = tax_level_list   

taxonomy_level_column_creation(taxsplit_otus, 0, "kingdom")
taxonomy_level_column_creation(taxsplit_otus, 1, "phylum")
taxonomy_level_column_creation(taxsplit_otus, 2, "class")
taxonomy_level_column_creation(taxsplit_otus, 3, "order")
taxonomy_level_column_creation(taxsplit_otus, 4, "family")
taxonomy_level_column_creation(taxsplit_otus, 5, "genus")
taxonomy_level_column_creation(taxsplit_otus, 6, "species")

#NEW as of 03.09.20: removing unclassified sections (e.g. s__)
#yes, this looks weird. Just removing "s__" didn't work, for unknown reasons. This is more generalizable anyway
def remove_blanks(otu_table, colname):
    for elem in otu_table[str(colname)]:
        if len(elem) <=4:
            dummy = elem
    otu_table[str(colname)].replace({dummy: ""}, inplace = True)
    del dummy
    del elem
    return

#Ahem, today I wrote a function without googling ANYTHING, thank you very much.

#remove_blanks(taxsplit_otus, "kingdom") can't actually do this one b/c of "UnboundLocalError: local variable 'dummy' referenced before assignment". Don't need it anyway.
remove_blanks(taxsplit_otus, "phylum")
remove_blanks(taxsplit_otus, "class")
remove_blanks(taxsplit_otus, "order")
remove_blanks(taxsplit_otus, "family")
remove_blanks(taxsplit_otus, "genus")
remove_blanks(taxsplit_otus, "species")    

#yes yes, very good. Now, to create columns where I stitch these back together. Preferrably in function form.
#with the 03.09.20 addition of the blanked cells, I need to figure out some sort of conditional, where stuff is only added if 
#the cell isn't blank.
#basically, if taxsplit_otus (next) is blank, don't add anything. Otherwise, do the command
#this is hacky but it works (I think?)
#kind of works. We still have the issue of multiple semi-colons

taxsplit_otus["phylum_level_sort"] = np.where(len(taxsplit_otus["phylum"]) > 1, (taxsplit_otus["kingdom"] + ";" + taxsplit_otus['phylum']), "")
taxsplit_otus["class_level_sort"] = np.where(len(taxsplit_otus["class"]) > 1, (taxsplit_otus["phylum_level_sort"] + ";" + taxsplit_otus['class']), "")
taxsplit_otus["order_level_sort"] = np.where(len(taxsplit_otus["order"]) > 1, (taxsplit_otus["class_level_sort"] + ";" + taxsplit_otus['order']), "")
taxsplit_otus["family_level_sort"] = np.where(len(taxsplit_otus["family"]) > 1, (taxsplit_otus["order_level_sort"] + ";" + taxsplit_otus['family']), "")
taxsplit_otus["genus_level_sort"] = np.where(len(taxsplit_otus["genus"]) > 1, (taxsplit_otus["family_level_sort"] + ";" + taxsplit_otus['genus']), "")
taxsplit_otus["species_level_sort"] = np.where(len(taxsplit_otus["species"]) > 1, (taxsplit_otus["genus_level_sort"] + ";" + taxsplit_otus['species']), "")

#fixed it!
taxsplit_otus["phylum_level_sort"] = taxsplit_otus["phylum_level_sort"].str.strip(' ;')
taxsplit_otus["class_level_sort"] = taxsplit_otus["class_level_sort"].str.strip(' ;')
taxsplit_otus["order_level_sort"] = taxsplit_otus["order_level_sort"].str.strip(' ;')
taxsplit_otus["family_level_sort"] = taxsplit_otus["family_level_sort"].str.strip(' ;')
taxsplit_otus["genus_level_sort"] = taxsplit_otus["genus_level_sort"].str.strip(' ;')
taxsplit_otus["species_level_sort"] = taxsplit_otus["species_level_sort"].str.strip(' ;')

#Added 04.15.20
taxsplit_otus["phylum_level_sort"] = taxsplit_otus["phylum_level_sort"].str.strip()
taxsplit_otus["class_level_sort"] = taxsplit_otus["class_level_sort"].str.strip()
taxsplit_otus["order_level_sort"] = taxsplit_otus["order_level_sort"].str.strip()
taxsplit_otus["family_level_sort"] = taxsplit_otus["family_level_sort"].str.strip()
taxsplit_otus["genus_level_sort"] = taxsplit_otus["genus_level_sort"].str.strip()
taxsplit_otus["species_level_sort"] = taxsplit_otus["species_level_sort"].str.strip()

#Now I have an OTU table that I can use for sorting, filtering, etc. I will use this moving forward

############################Step 1: remove all obvious contaminants############################
#here is my contaminants list:
contaminants = pd.read_csv("06.16.20_Burns_contaminant_list.txt", sep="    ", header = None, engine='python') #new contaminants list
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
        
        
#then, create a bunch of indices where the contaminants list is cross-referenced with taxsplit_otus. This will basically create a reference index to pull out
#contaminant taxa from taxsplit_otus
kingdom_index = taxsplit_otus[taxsplit_otus.kingdom.isin(kingdom_contaminants)].index  #pulls out indices (here OTUs) where the
#corresponding value in the kingdom column is in the list kingdom_contaminants
taxsplit_otus.drop(kingdom_index, inplace=True)
        #repeat for all levels
phylum_index = taxsplit_otus[taxsplit_otus.phylum_level_sort.isin(phylum_contaminants)].index  
taxsplit_otus.drop(phylum_index, inplace=True)

class_index = taxsplit_otus[taxsplit_otus.class_level_sort.isin(class_contaminants)].index  
taxsplit_otus.drop(class_index, inplace=True)

order_index = taxsplit_otus[taxsplit_otus.order_level_sort.isin(order_contaminants)].index  
taxsplit_otus.drop(order_index, inplace=True)

family_index = taxsplit_otus[taxsplit_otus.family_level_sort.isin(family_contaminants)].index  
taxsplit_otus.drop(family_index, inplace=True)

genus_index = taxsplit_otus[taxsplit_otus.genus_level_sort.isin(genus_contaminants)].index  
taxsplit_otus.drop(genus_index, inplace=True)

#TA-DA!!!
#save the filtered dataframe (full) to csv
#taxsplit_otus.to_csv("04.27.20_Burns_2015_OTUS_contaminants_removed.csv")
#this is the contaminants removed dataset. No abundance filtering has been performed here. 

#####New as of 04.15.20: relative abundance filtering on whole OTU table########

#step 0: pull all columns that need to be summed (i.e. all OTU columns)- will use a lot
sample_names = []            #this is a list of all the column names that will need to be scanned for values >0.001
for col in taxsplit_otus.columns:       #shouldn't matter which dataframe we use here but use raw one just in case
    if 'Sample_' in str(col):
        sample_names.append(str(col))
    del col

#step 1: collapse by class- not really applicable here but the code won't run without it
otu_table_collapsed = taxsplit_otus.groupby(['taxonomy_orig', '#OTU_ID'])[sample_names].sum()
#fine, but we actually need everything else in taxsplit_otus that's not in sample_names
#or IS THERE???
#maybe we just keep taxonomy_orig and #OTU_ID. the rest of the taxonomy stuff we can add again later.

#step 2: calculate relative abundances
otu_table_relabund = otu_table_collapsed .loc[:,].div(otu_table_collapsed .sum(axis=0))
#it is worth noting that this works BECAUSE in the groupby function in step 1, taxonomy_orig and #OTU_ID are both
#set as the index. Otherwise, can't do maths on strings (which is the "ValueError: operands could not be broadcast 
#together with shapes (1192,) (8344, )" error I was getting for a while

#step 4: filter by abundance (at least 0.0001 in at least 20% of samples)
#define the number of samples that hit the cutoff for each taxon
num_above_cutoff = [(otu_table_relabund[sample_names] >= 0.001).sum(1)]

#convert to a list
for elem in num_above_cutoff:
    otu_table_num_above_cutoff = elem.values.tolist()       #this converts a series to a list
    
#add this list as a column
otu_table_relabund['Num_above_cutoff'] = otu_table_num_above_cutoff

#step 5: calculate the number of samples in the total dataframe- will need this value later
numsamples =0
for col in otu_table_relabund:
    if 'Sample_' in str(col):
        numsamples = numsamples+1
        
#step 6: apply the filter
otu_table_relabund2 = otu_table_relabund.loc[otu_table_relabund['Num_above_cutoff'] >= (numsamples*0.2)]
#also do filtering at 10% of samples
otu_table_relabund3 = otu_table_relabund.loc[otu_table_relabund['Num_above_cutoff'] >= (numsamples*0.1)]
#NEW as of 06.16.20: find out how many samples were abundance-filtered out. 
otu_table_relabund4 = otu_table_relabund.loc[otu_table_relabund['Num_above_cutoff'] < (numsamples*0.1)]

#step 7: reset the indices
otu_table_relabund2.reset_index(inplace=True)
otu_table_relabund3.reset_index(inplace=True)
#NEW as of 06.16.20: find out how many samples were abundance-filtered out. 
otu_table_relabund4.reset_index(inplace=True)

#step 8: redo all the phylogeny generation stuff from above. 
#have to copy the stupid dataframes or get 'settingwithcopywarning' errors
Burns_otu_table_20 = copy.deepcopy(otu_table_relabund2)
Burns_otu_table_10 = copy.deepcopy(otu_table_relabund3)
#NEW as of 06.16.20: find out how many samples were abundance-filtered out. 
Burns_otu_table_lowab = copy.deepcopy(otu_table_relabund4)

Burns_otu_table_20['taxonomy'] = Burns_otu_table_20['taxonomy_orig']     #copy the taxonomy column to save it for later
Burns_otu_table_20["taxonomy"] = Burns_otu_table_20["taxonomy"].str.split(";", expand=False)    #split the values in the taxonomy column into a list; delimiter is a semicolon

Burns_otu_table_10['taxonomy'] = Burns_otu_table_10['taxonomy_orig']     #copy the taxonomy column to save it for later
Burns_otu_table_10["taxonomy"] = Burns_otu_table_10["taxonomy"].str.split(";", expand=False)    #split the values in the taxonomy column into a list; delimiter is a semicolon

#NEW as of 06.16.20: find out how many samples were abundance-filtered out. 
Burns_otu_table_lowab['taxonomy'] = Burns_otu_table_lowab['taxonomy_orig']     #copy the taxonomy column to save it for later
Burns_otu_table_lowab["taxonomy"] = Burns_otu_table_lowab["taxonomy"].str.split(";", expand=False)    #split the values in the taxonomy column into a list; delimiter is a semicolon

#now, to create new columns based on different taxonomic levels
#this function was defined above 
taxonomy_level_column_creation(Burns_otu_table_20, 0, "kingdom")
taxonomy_level_column_creation(Burns_otu_table_20, 1, "phylum")
taxonomy_level_column_creation(Burns_otu_table_20, 2, "class")
taxonomy_level_column_creation(Burns_otu_table_20, 3, "order")
taxonomy_level_column_creation(Burns_otu_table_20, 4, "family")
taxonomy_level_column_creation(Burns_otu_table_20, 5, "genus")
taxonomy_level_column_creation(Burns_otu_table_20, 6, "species")

taxonomy_level_column_creation(Burns_otu_table_10, 0, "kingdom")
taxonomy_level_column_creation(Burns_otu_table_10, 1, "phylum")
taxonomy_level_column_creation(Burns_otu_table_10, 2, "class")
taxonomy_level_column_creation(Burns_otu_table_10, 3, "order")
taxonomy_level_column_creation(Burns_otu_table_10, 4, "family")
taxonomy_level_column_creation(Burns_otu_table_10, 5, "genus")
taxonomy_level_column_creation(Burns_otu_table_10, 6, "species")

#NEW as of 06.16.20: find out how many samples were abundance-filtered out. 
taxonomy_level_column_creation(Burns_otu_table_lowab, 0, "kingdom")
taxonomy_level_column_creation(Burns_otu_table_lowab, 1, "phylum")
taxonomy_level_column_creation(Burns_otu_table_lowab, 2, "class")
taxonomy_level_column_creation(Burns_otu_table_lowab, 3, "order")
taxonomy_level_column_creation(Burns_otu_table_lowab, 4, "family")
taxonomy_level_column_creation(Burns_otu_table_lowab, 5, "genus")
taxonomy_level_column_creation(Burns_otu_table_lowab, 6, "species")


#NEW as of 03.09.20: removing unclassified sections (e.g. s__)
#remove_blanks function is defined above

#yes yes, very good. Now, to create columns where I stitch these back together. Preferrably in function form.
#with the 03.09.20 addition of the blanked cells, I need to figure out some sort of conditional, where stuff is only added if 
#the cell isn't blank.
#basically, if taxsplit_otus (next) is blank, don't add anything. Otherwise, do the command
#this is hacky but it works (I think?)
#kind of works. We still have the issue of multiple semi-colons

Burns_otu_table_20["phylum_level_sort"] = np.where(len(Burns_otu_table_20["phylum"]) > 1, (Burns_otu_table_20["kingdom"] + ";" + Burns_otu_table_20['phylum']), "")
Burns_otu_table_20["class_level_sort"] = np.where(len(Burns_otu_table_20["class"]) > 1, (Burns_otu_table_20["phylum_level_sort"] + ";" + Burns_otu_table_20['class']), "")
Burns_otu_table_20["order_level_sort"] = np.where(len(Burns_otu_table_20["order"]) > 1, (Burns_otu_table_20["class_level_sort"] + ";" + Burns_otu_table_20['order']), "")
Burns_otu_table_20["family_level_sort"] = np.where(len(Burns_otu_table_20["family"]) > 1, (Burns_otu_table_20["order_level_sort"] + ";" + Burns_otu_table_20['family']), "")
Burns_otu_table_20["genus_level_sort"] = np.where(len(Burns_otu_table_20["genus"]) > 1, (Burns_otu_table_20["family_level_sort"] + ";" + Burns_otu_table_20['genus']), "")
Burns_otu_table_20["species_level_sort"] = np.where(len(Burns_otu_table_20["species"]) > 1, (Burns_otu_table_20["genus_level_sort"] + ";" + Burns_otu_table_20['species']), "")

#fixed it!
Burns_otu_table_20["phylum_level_sort"] = Burns_otu_table_20["phylum_level_sort"].str.strip(' ;')
Burns_otu_table_20["class_level_sort"] = Burns_otu_table_20["class_level_sort"].str.strip(' ;')
Burns_otu_table_20["order_level_sort"] = Burns_otu_table_20["order_level_sort"].str.strip(' ;')
Burns_otu_table_20["family_level_sort"] = Burns_otu_table_20["family_level_sort"].str.strip(' ;')
Burns_otu_table_20["genus_level_sort"] = Burns_otu_table_20["genus_level_sort"].str.strip(' ;')
Burns_otu_table_20["species_level_sort"] = Burns_otu_table_20["species_level_sort"].str.strip(' ;')

#Added 04.15.20
Burns_otu_table_20["phylum_level_sort"] = Burns_otu_table_20["phylum_level_sort"].str.strip()
Burns_otu_table_20["class_level_sort"] = Burns_otu_table_20["class_level_sort"].str.strip()
Burns_otu_table_20["order_level_sort"] = Burns_otu_table_20["order_level_sort"].str.strip()
Burns_otu_table_20["family_level_sort"] = Burns_otu_table_20["family_level_sort"].str.strip()
Burns_otu_table_20["genus_level_sort"] = Burns_otu_table_20["genus_level_sort"].str.strip()
Burns_otu_table_20["species_level_sort"] = Burns_otu_table_20["species_level_sort"].str.strip()

Burns_otu_table_10["phylum_level_sort"] = np.where(len(Burns_otu_table_10["phylum"]) > 1, (Burns_otu_table_10["kingdom"] + ";" + Burns_otu_table_10['phylum']), "")
Burns_otu_table_10["class_level_sort"] = np.where(len(Burns_otu_table_10["class"]) > 1, (Burns_otu_table_10["phylum_level_sort"] + ";" + Burns_otu_table_10['class']), "")
Burns_otu_table_10["order_level_sort"] = np.where(len(Burns_otu_table_10["order"]) > 1, (Burns_otu_table_10["class_level_sort"] + ";" + Burns_otu_table_10['order']), "")
Burns_otu_table_10["family_level_sort"] = np.where(len(Burns_otu_table_10["family"]) > 1, (Burns_otu_table_10["order_level_sort"] + ";" + Burns_otu_table_10['family']), "")
Burns_otu_table_10["genus_level_sort"] = np.where(len(Burns_otu_table_10["genus"]) > 1, (Burns_otu_table_10["family_level_sort"] + ";" + Burns_otu_table_10['genus']), "")
Burns_otu_table_10["species_level_sort"] = np.where(len(Burns_otu_table_10["species"]) > 1, (Burns_otu_table_10["genus_level_sort"] + ";" + Burns_otu_table_10['species']), "")

#fixed it!
Burns_otu_table_10["phylum_level_sort"] = Burns_otu_table_10["phylum_level_sort"].str.strip(' ;')
Burns_otu_table_10["class_level_sort"] = Burns_otu_table_10["class_level_sort"].str.strip(' ;')
Burns_otu_table_10["order_level_sort"] = Burns_otu_table_10["order_level_sort"].str.strip(' ;')
Burns_otu_table_10["family_level_sort"] = Burns_otu_table_10["family_level_sort"].str.strip(' ;')
Burns_otu_table_10["genus_level_sort"] = Burns_otu_table_10["genus_level_sort"].str.strip(' ;')
Burns_otu_table_10["species_level_sort"] = Burns_otu_table_10["species_level_sort"].str.strip(' ;')

#Added 04.15.20
Burns_otu_table_10["phylum_level_sort"] = Burns_otu_table_10["phylum_level_sort"].str.strip()
Burns_otu_table_10["class_level_sort"] = Burns_otu_table_10["class_level_sort"].str.strip()
Burns_otu_table_10["order_level_sort"] = Burns_otu_table_10["order_level_sort"].str.strip()
Burns_otu_table_10["family_level_sort"] = Burns_otu_table_10["family_level_sort"].str.strip()
Burns_otu_table_10["genus_level_sort"] = Burns_otu_table_10["genus_level_sort"].str.strip()
Burns_otu_table_10["species_level_sort"] = Burns_otu_table_10["species_level_sort"].str.strip()

#NEW as of 06.16.20: find out how many samples were abundance-filtered out. 
Burns_otu_table_lowab["phylum_level_sort"] = np.where(len(Burns_otu_table_lowab["phylum"]) > 1, (Burns_otu_table_lowab["kingdom"] + ";" + Burns_otu_table_lowab['phylum']), "")
Burns_otu_table_lowab["class_level_sort"] = np.where(len(Burns_otu_table_lowab["class"]) > 1, (Burns_otu_table_lowab["phylum_level_sort"] + ";" + Burns_otu_table_lowab['class']), "")
Burns_otu_table_lowab["order_level_sort"] = np.where(len(Burns_otu_table_lowab["order"]) > 1, (Burns_otu_table_lowab["class_level_sort"] + ";" + Burns_otu_table_lowab['order']), "")
Burns_otu_table_lowab["family_level_sort"] = np.where(len(Burns_otu_table_lowab["family"]) > 1, (Burns_otu_table_lowab["order_level_sort"] + ";" + Burns_otu_table_lowab['family']), "")
Burns_otu_table_lowab["genus_level_sort"] = np.where(len(Burns_otu_table_lowab["genus"]) > 1, (Burns_otu_table_lowab["family_level_sort"] + ";" + Burns_otu_table_lowab['genus']), "")
Burns_otu_table_lowab["species_level_sort"] = np.where(len(Burns_otu_table_lowab["species"]) > 1, (Burns_otu_table_lowab["genus_level_sort"] + ";" + Burns_otu_table_lowab['species']), "")

Burns_otu_table_lowab["phylum_level_sort"] = Burns_otu_table_lowab["phylum_level_sort"].str.strip(' ;')
Burns_otu_table_lowab["class_level_sort"] = Burns_otu_table_lowab["class_level_sort"].str.strip(' ;')
Burns_otu_table_lowab["order_level_sort"] = Burns_otu_table_lowab["order_level_sort"].str.strip(' ;')
Burns_otu_table_lowab["family_level_sort"] = Burns_otu_table_lowab["family_level_sort"].str.strip(' ;')
Burns_otu_table_lowab["genus_level_sort"] = Burns_otu_table_lowab["genus_level_sort"].str.strip(' ;')
Burns_otu_table_lowab["species_level_sort"] = Burns_otu_table_lowab["species_level_sort"].str.strip(' ;')

Burns_otu_table_lowab["phylum_level_sort"] = Burns_otu_table_lowab["phylum_level_sort"].str.strip()
Burns_otu_table_lowab["class_level_sort"] = Burns_otu_table_lowab["class_level_sort"].str.strip()
Burns_otu_table_lowab["order_level_sort"] = Burns_otu_table_lowab["order_level_sort"].str.strip()
Burns_otu_table_lowab["family_level_sort"] = Burns_otu_table_lowab["family_level_sort"].str.strip()
Burns_otu_table_lowab["genus_level_sort"] = Burns_otu_table_lowab["genus_level_sort"].str.strip()
Burns_otu_table_lowab["species_level_sort"] = Burns_otu_table_lowab["species_level_sort"].str.strip()


Burns_otu_table_20.to_csv("06.16.20_otu_table_abundance_filtered_0.2.csv")    
Burns_otu_table_10.to_csv("06.16.20_otu_table_abundance_filtered_0.1.csv")

Burns_otu_table_lowab.to_csv("06.16.20_otu_table_below_0.1_abundance_taxa.csv")

        
############################################################################################################################################
#########################################Relative abundance filtering by taxonomic level########################################################
############################################################################################################################################


############################Step 2: create separate dataframes to filter/collapse at each taxonomic level############################
#get a list of the column names and their indices- will need these (probably) to make new dataframes
col_names_list = []
for col in taxsplit_otus.columns:
    col_names_list.append(str(col))
    del col
#interesting, I ended up with 100 columns. 0 is the OTU_ID, 1-88 inclusive are the samples,
#89 is the semicolon-seperated taxonomy list, 90 is the original taxonomy string, 91-97 inclusive are the
#individual taxonomic levels, 98-102 are the rebound taxonomies at the phylum to genus level
    
#create seperate dataframes 
kingdom_level_otus = taxsplit_otus.drop(["taxonomy",  v"phylum", "class", "order", "family", "genus", "species", "phylum_level_sort", 
                                         "class_level_sort", "order_level_sort", "family_level_sort", "genus_level_sort"], axis = 'columns')    
phylum_level_otus = taxsplit_otus.drop(["taxonomy", "kingdom", "class", "order", "family", "genus", "species", 
                                        "class_level_sort", "order_level_sort", "family_level_sort", "genus_level_sort"], axis = 'columns')    
class_level_otus = taxsplit_otus.drop(["taxonomy", "kingdom", "phylum", "order", "family", "genus", "species", "phylum_level_sort", 
                                        "order_level_sort", "family_level_sort", "genus_level_sort"], axis = 'columns')    
order_level_otus = taxsplit_otus.drop(["taxonomy", "kingdom", "phylum", "class", "family", "genus", "species", "phylum_level_sort", 
                                       "class_level_sort", "family_level_sort", "genus_level_sort"], axis = 'columns')    
family_level_otus = taxsplit_otus.drop(["taxonomy", "kingdom", "phylum", "class", "order", "genus", "species", "phylum_level_sort", 
                                        "class_level_sort", "order_level_sort", "genus_level_sort"], axis = 'columns')    
genus_level_otus = taxsplit_otus.drop(["taxonomy", "kingdom", "phylum", "class", "order", "family", "species", "phylum_level_sort", 
                                       "class_level_sort", "order_level_sort", "family_level_sort"], axis = 'columns')    

############################Step 3: collapse at each taxonomic level############################

#note that kingdom level sorting only involves removing archaea- that was done in the contaminant removal step above. So kingdom code is short
kingdom_level_otus.sort_values(by=["kingdom"])

#step 0: pull all columns that need to be summed (i.e. all OTU columns)- will use a lot
sample_names = []            #this is a list of all the column names that will need to be scanned for values >0.001
for col in taxsplit_otus.columns:       #shouldn't matter which dataframe we use here but use raw one just in case
    if 'Sample_' in str(col):
        sample_names.append(str(col))

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
    tax_level_relabund2.to_csv("_relabund_0.1.csv")  
    for file in os.listdir():                       #rename files by their appropriate diet; couldn't figure out a way to do this in the above code
        src=file
        if fnmatch.fnmatch(file, "_relabund_0.1.csv"):
            dst = str(tax_level)+str('_level')+file
            os.rename(src,dst)
    return

tax_level_filtering("phylum", phylum_level_otus, 'phylum_level_sort')
tax_level_filtering("class", class_level_otus, 'class_level_sort')
tax_level_filtering("order", order_level_otus, 'order_level_sort')
tax_level_filtering("family", family_level_otus, 'family_level_sort')
tax_level_filtering("genus", genus_level_otus, 'genus_level_sort')




#############JUNK
#have a conditional way to do this now
taxsplit_otus["phylum_level_sort"] = taxsplit_otus['kingdom'] + "; " + taxsplit_otus['phylum'] #this one needs to be done first
taxsplit_otus["class_level_sort"] = taxsplit_otus["phylum_level_sort"] + "; " + taxsplit_otus['class']
taxsplit_otus["order_level_sort"] = taxsplit_otus["class_level_sort"] + "; " + taxsplit_otus['order']
taxsplit_otus["family_level_sort"] = taxsplit_otus["order_level_sort"] + "; " + taxsplit_otus['family']
taxsplit_otus["genus_level_sort"] = taxsplit_otus["family_level_sort"] + "; " + taxsplit_otus['genus']
#I could probably do this iteratively but right now I am unwell and can't think of how to do it. The code
#should work for any OTU table anyway, especially if I name the dataframe taxsplit_otus.













