#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 16:29:25 2021

@author: adamo010
"""

import scipy
import numpy as np
import csv
import subprocess
import os
import pandas as pd
import fnmatch
import shutil
import glob
import re
import copy

#you know, I don't actually have a good taxonomy key
preprocessed_taxonomy = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal/Burns_absolute_OTU_table_phylogeny_corrected.csv")
preprocessed_taxonomy.rename(columns={"#OTU_ID":"OTU_ID"}, inplace = True)
preprocessed_taxonomy.drop(preprocessed_taxonomy.loc[:, 'Sample_27':'Sample_42'].columns, axis = 1, inplace=True) 
preprocessed_taxonomy.to_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Burns_2015_OTU_id_taxonomy_table.csv")

#ugh, how do I expand this again?
taxsplit_otus = copy.deepcopy(preprocessed_taxonomy)   #copy for safety
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

taxsplit_otus.drop("taxonomy", axis=1, inplace=True)

#remove "k__", "g__", etc. from each column. This is the first three characters for kingdom and the first four for all otehrs, so just do all those
taxsplit_otus["kingdom"] = taxsplit_otus["kingdom"].str[3:]
taxsplit_otus["phylum"] = taxsplit_otus["phylum"].str[4:]
taxsplit_otus["class"] = taxsplit_otus["class"].str[4:]
taxsplit_otus["order"] = taxsplit_otus["order"].str[4:]
taxsplit_otus["family"] = taxsplit_otus["family"].str[4:]
taxsplit_otus["genus"] = taxsplit_otus["genus"].str[4:]
taxsplit_otus["species"] = taxsplit_otus["species"].str[4:]

taxsplit_otus.to_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Burns_2015_OTU_id_taxonomy_table_tax_levels_separated.csv")


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
