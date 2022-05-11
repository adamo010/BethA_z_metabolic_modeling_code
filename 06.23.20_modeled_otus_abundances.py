#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 16:39:06 2020

@author: adamo010
"""

#now, what I would like to do is some abundance calculations

#need to create an absolute abundance OTU table from just the OTUs that have models

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
import ast

otus_with_models = pd.read_csv("06.23.20_summary_of_otus_with_models.csv")
preprocessed_otus = pd.read_csv("Contaminant_removal/Burns_absolute_OTU_table.csv")

#STEP 1: borrow from 06.16.20_OTU_table_filtering_pipeline to get OTU table up to snuff

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

#removing unclassified sections (e.g. s__)
#yes, this looks weird. Just removing "s__" didn't work, for unknown reasons. This is more generalizable anyway
def remove_blanks(otu_table, colname):
    for elem in otu_table[str(colname)]:
        if len(elem) <=4:
            dummy = elem
    otu_table[str(colname)].replace({dummy: ""}, inplace = True)
    del dummy
    del elem
    return

#remove_blanks(taxsplit_otus, "kingdom") can't actually do this one b/c of "UnboundLocalError: local variable 'dummy' referenced before assignment". Don't need it anyway.
remove_blanks(taxsplit_otus, "phylum")
remove_blanks(taxsplit_otus, "class")
remove_blanks(taxsplit_otus, "order")
remove_blanks(taxsplit_otus, "family")
remove_blanks(taxsplit_otus, "genus")
remove_blanks(taxsplit_otus, "species")    

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

#STEP 2. Filter taxsplit OTUs based on otus present in otus_with_models
#get a list from the OTU_ID column in otus_with_models
otus_with_models_list = otus_with_models['OTU_ID'].tolist()

#now filter the OTU table based on that list
filtered_otu_table = taxsplit_otus[taxsplit_otus['#OTU_ID'].isin(otus_with_models_list)]

filtered_otu_table.to_csv("06.23.20_Burns_absolute_abundances_filtered_by_has_model.csv")
#be nice to have this guy kicking around for a bit

#######STEP 3: calculate relative abundances

#step 0: pull all columns that need to be summed (i.e. all OTU columns)- will use a lot
sample_names = []            #this is a list of all the column names that will need to be scanned for values >0.001
for col in filtered_otu_table.columns:       #shouldn't matter which dataframe we use here but use raw one just in case
    if 'Sample_' in str(col):
        sample_names.append(str(col))
    del col

#step 1: collapse by class- not really applicable here but the code won't run without it
otu_table_collapsed = filtered_otu_table.groupby(['taxonomy_orig', '#OTU_ID'])[sample_names].sum()
#fine, but we actually need everything else in taxsplit_otus that's not in sample_names
#or IS THERE???
#maybe we just keep taxonomy_orig and #OTU_ID. the rest of the taxonomy stuff we can add again later.

#step 2: calculate relative abundances
otu_table_relabund = otu_table_collapsed .loc[:,].div(otu_table_collapsed .sum(axis=0))
#it is worth noting that this works BECAUSE in the groupby function in step 1, taxonomy_orig and #OTU_ID are both
#set as the index. Otherwise, can't do maths on strings (which is the "ValueError: operands could not be broadcast 
#together with shapes (1192,) (8344, )" error I was getting for a while

#well... okay. Now need to think about how to incorporate the CRC vs. healthy tissues in here. Are there abundance differences?
#so, basically, need to do some real microbiome work. 

















