#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 11:53:54 2021

@author: adamo010
"""

#for Niccolai2020 dataset

#based on 07.20.21_Picking_AGORA_and_MAMBO_models_single_tax_levels_Hale2018_data.py
#and 07.29.21_Picking_Carved_models_Hale2018_data.py

#basically, the goal is to taxonomically match all OTUs for whom I have now (theoretically) Carved models
#so I can have 100% matches for the 49 OTUs which lack MAMBO or AGORA models.

#using a new database today- the CarveMe database. Dont remember how I did this last time, but for Hale2018
#data this is how it's going to go.

import scipy
import numpy as np
import cobra
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

#step 1: import data.
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Niccolai_2020_model_picking")

carveme = pd.read_csv("10.19.21_Niccolai2020_Carved_model_taxonomy.csv")
otus = pd.read_csv("10.19.21_Niccolai2020_OTUs_without_mambo_models_single_tax_level_hand_picks_V3.csv") 

#really need to clear out excess columns- there are too many
sorted(otus) #prints column names
otus_clean = otus.drop(columns=["Unnamed: 0", 'agora_family_matched_model_list', 'agora_genus_matched_model_list',
       'agora_genus2_matched_model_list', 'all_agora_matched_models',
       'mambo_family_matched_model_list', 'mambo_genus_new_matched_model_list',
       'mambo_genus_extracted_matched_model_list', 'all_mambo_matched_models'])
#columns of interest: Kingdom, Phylum, Class, Order, Family, Genus, Species
#columns in carveme:  Kingdom, Phylum, Class, Order, Family, Genus, Species
    
###############step 2a: create dictionaries for each taxonomic level carveme
    #columns in carveme: family, genus, genus2, species
#carveme family_
unique_carveme_family = carveme.Family.unique()
unique_carveme_family = pd.DataFrame(unique_carveme_family)    #import to pandas
unique_carveme_family['model_list'] = np.empty((len(unique_carveme_family), 0)).tolist()  #add an empty 'model list' column
unique_carveme_family.reset_index()
unique_carveme_family= unique_carveme_family.rename(columns={0: "family_sort"})       #rename column
carveme_family_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_carveme_family['family_sort']:  
    family_name_key = str(row)
    temp_df = carveme[carveme["Family"].str.contains(str(row))]        #temp_df should now contain only carveme rows that match that row in unique_carveme_genera2 
    templist = []
    for row in temp_df['model_file_name']:
        templist.append(str(row))
    del temp_df
    carveme_family_models_dict[family_name_key] = templist
del unique_carveme_family, family_name_key, row, templist
#carveme genus    
unique_carveme_genus = carveme.Genus.unique()
unique_carveme_genus = pd.DataFrame(unique_carveme_genus)    #import to pandas
unique_carveme_genus['model_list'] = np.empty((len(unique_carveme_genus), 0)).tolist()  #add an empty 'model list' column
unique_carveme_genus.reset_index()
unique_carveme_genus= unique_carveme_genus.rename(columns={0: "genus_sort"})       #rename column
carveme_genus_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_carveme_genus['genus_sort']:  
    genus_name_key = str(row)
    temp_df = carveme[carveme["Genus"].str.contains(str(row))]        #temp_df should now contain only carveme rows that match that row in unique_carveme_genera2 
    templist = []
    for row in temp_df['model_file_name']:
        templist.append(str(row))
    del temp_df
    carveme_genus_models_dict[genus_name_key] = templist
del unique_carveme_genus, genus_name_key, row, templist
#carveme species: NOT NEEDED HERE
#unique_carveme_species = carveme.Species.unique()
#unique_carveme_species = pd.DataFrame(unique_carveme_species)    #import to pandas
#unique_carveme_species['model_list'] = np.empty((len(unique_carveme_species), 0)).tolist()  #add an empty 'model list' column
#unique_carveme_species.reset_index()
#unique_carveme_species= unique_carveme_species.rename(columns={0: "species_sort"})       #rename column
#carveme_species_models_dict ={}         #this is where I'll store the key-value pairs for species- [list of models]
#for row in unique_carveme_species['species_sort']:  
    #species_name_key = str(row)
    #temp_df = carveme[carveme["Species"].str.contains(str(row))]        #temp_df should now contain only carveme rows that match that row in unique_carveme_genera2 
    #templist = []
    #for row in temp_df['model_file_name']:
        #templist.append(str(row))
    #del temp_df
    #carveme_species_models_dict[species_name_key] = templist
#del unique_carveme_species, species_name_key, row, templist

###############step 3: create new otu tables based on 'id to species level', 'id to genus level', 'id to family level' without replacement
#note that this section is generalized; don't need to repeat for MAMBO.
# "fgs_level_sort", "fg_level_sort", "gs_level_sort", "fgs2_level_sort", "fg2_level_sort", "gs2_level_sort"
otus_plus_tax = copy.deepcopy(otus_clean)
#now bin by tax level
#otus_to_species_level = otus_plus_tax[otus_plus_tax['Species'].notnull()]   #hey, look at that
otus_to_genus_level = otus_plus_tax[otus_plus_tax['Genus'].notnull()]   #you may award me the genuis medal now
otus_to_family_level = otus_plus_tax[otus_plus_tax['Family'].notnull() & otus_plus_tax['Genus'].isnull()]   
#all right, now things are gonna get a bit weird
#can I get away with iterating through otus_plus_tax and just... adding whatever models fit? And keeping the col with the biggest # models?

#############step 4a: carveme: map models and reorganize columns
otus_plus_tax_plus_models = copy.deepcopy(otus_plus_tax)
otus_plus_tax_plus_models["carveme_family_matched_model_list"] = otus_plus_tax_plus_models["Family"].map(carveme_family_models_dict)
otus_plus_tax_plus_models["carveme_genus_matched_model_list"] = otus_plus_tax_plus_models["Genus"].map(carveme_genus_models_dict)
#otus_plus_tax_plus_models["carveme_species_matched_model_list"] = otus_plus_tax_plus_models["Species"].map(carveme_species_models_dict)    

#first, fill in nan values with blank lists.
#make a list of the columns which need nans replaced with blank lists
list_of_cols = ["carveme_family_matched_model_list", "carveme_genus_matched_model_list"]

#for each of these olumns, replace all nans with [] (blank lists)
for col in list_of_cols:
    otus_plus_tax_plus_models.loc[otus_plus_tax_plus_models[col].isnull(),[col]] = \
    otus_plus_tax_plus_models.loc[otus_plus_tax_plus_models[col].isnull(),col].apply(lambda x: [])

#create a new column. all_carveme_matched_models, that contains a merged list of all matched models at all levels
otus_plus_tax_plus_models["all_carveme_matched_models"] = otus_plus_tax_plus_models["carveme_family_matched_model_list"] + \
    otus_plus_tax_plus_models["carveme_genus_matched_model_list"]

#this addition created some duplicates (e.g. where gs and gs2 pulled the same models)    
otus_plus_tax_plus_models['all_carveme_matched_models'] = otus_plus_tax_plus_models['all_carveme_matched_models'].apply(lambda x: list(set(x)))

#nice. Now, we move onto MAMBO. But first, split this datframe into "has carveme model" and "needs mambo model"
#note that we can't use the 'null' thing because there are blank lists in this column
otus_with_carveme_models = otus_plus_tax_plus_models[~otus_plus_tax_plus_models['all_carveme_matched_models'].str.len().eq(0)]  
otus_without_carveme_models = otus_plus_tax_plus_models[otus_plus_tax_plus_models['all_carveme_matched_models'].str.len().eq(0)]  
#NOTE! otus_without_carveme_models SHOULD be empty, as we created the models specifically to match the remaining OTUs

#SAVE
otus_with_carveme_models.to_csv("10.19.21_Niccolai2020_OTUs_with_carveme_models.csv")


