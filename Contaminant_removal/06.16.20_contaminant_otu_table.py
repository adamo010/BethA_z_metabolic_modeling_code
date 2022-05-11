#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 17:59:22 2020

@author: adamo010
"""
#update from 05.28.20: new contaminants list to remove duplicates

#so, apparently, when I filtered all the contaminants out, I didn't think to save those contaminants as a csv file.

#borrow heavily from 04.27.20_OTU_table_filtering_pipeline_clean
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

taxsplit_otus.rename(columns = {"#OTU_ID": "OTU_ID"}, inplace=True)

#find all indices that correspond to taxonomy level-specific contaminants
#filter the taxsplit OTU table to only include those indices

genus_index = taxsplit_otus[taxsplit_otus.genus_level_sort.isin(genus_contaminants)].index 
genus_table  = taxsplit_otus[taxsplit_otus.index.isin(genus_index)]

family_index = taxsplit_otus[taxsplit_otus.family_level_sort.isin(family_contaminants)].index  
family_table  = taxsplit_otus[taxsplit_otus.index.isin(family_index)]

order_index = taxsplit_otus[taxsplit_otus.order_level_sort.isin(order_contaminants)].index
order_table  = taxsplit_otus[taxsplit_otus.index.isin(order_index)]

class_index = taxsplit_otus[taxsplit_otus.class_level_sort.isin(class_contaminants)].index 
class_table  = taxsplit_otus[taxsplit_otus.index.isin(class_index)]

phylum_index = taxsplit_otus[taxsplit_otus.phylum_level_sort.isin(phylum_contaminants)].index
phylum_table  = taxsplit_otus[taxsplit_otus.index.isin(phylum_index)]

kingdom_index = taxsplit_otus[taxsplit_otus.kingdom.isin(kingdom_contaminants)].index 
kingdom_table  = taxsplit_otus[taxsplit_otus.index.isin(kingdom_index)]

frames =[genus_table, family_table, order_table, class_table, phylum_table, kingdom_table]
contaminant_otu_table = pd.concat(frames)

#TA-DA!!!
#save the filtered dataframe (full) to csv
contaminant_otu_table.to_csv("06.16.20_table_of_contaminant_OTUs.csv")
#this is the contaminants removed dataset. No abundance filtering has been performed here. 

