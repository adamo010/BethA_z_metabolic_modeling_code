#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 07:56:07 2021

@author: adamo010
"""
#The goal here is to match AGORA and MAMBO models to OTUs from Hale 2018 data.
#Specifically, I want to pick models with single taxonomic levels.
#e.g. pick based on genus only; this is to account for changes in taxonomies between
#when the AGORA/MAMBO model taxonomies were generated, and when Hale2018 OTU table was generated.

#all this is based on code from 07.19.21_Picking_AGORA_and_MAMBO_models_revised_taxonomies_Hale2018_data.py

import micom
import gurobipy
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
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Hale_2018_model_picking")

agora = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/AGORA_model_phylogeny_table.csv")
mambo = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/07.23.20_MAMBO_model_taxonomy_assigned_plus_Treponema.csv")
otus = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Hale_2018_model_picking/07.19.21_Hale2018_OTUs_without_mambo_or_agora_models_round2_picks_indiv_specific_comms.csv") #this is the new contaminant-free table

#really need to clear out excess columns- there are too many
sorted(otus)
otus_clean = otus.drop(columns=["Unnamed: 0", 'agora_fg2_matched_model_list', 'agora_fg_matched_model_list', 'agora_fgs2_matched_model_list',\
                                 'agora_fgs_matched_model_list', 'agora_gs2_matched_model_list', 'agora_gs_matched_model_list',\
                                     'all_agora_matched_models', 'all_mambo_matched_models', 'clean_species_level_sort',\
                                          'family_level_sort', 'fg_level_sort', 'fgs_level_sort', 'gs_level_sort', 'keep_or_toss',\
                                               'mambo_fg2_matched_model_list', 'mambo_fg_matched_model_list', 'mambo_fgs2_matched_model_list',\
                                                    'mambo_fgs_matched_model_list', 'mambo_gs2_matched_model_list', 'mambo_gs_matched_model_list',])
#columns of interest: Kingdom, Phylum, Class, Order, Family, Genus, Species
#columns in agora: family, genus, genus2, species
#columns in mambo: Family, Genus_new, Genus_extracted, Species_extracted
    
###############step 2a: create dictionaries for each taxonomic level AGORA
    #columns in agora: family, genus, genus2, species
#agora family_
unique_agora_family = agora.family.unique()
unique_agora_family = pd.DataFrame(unique_agora_family)    #import to pandas
unique_agora_family['model_list'] = np.empty((len(unique_agora_family), 0)).tolist()  #add an empty 'model list' column
unique_agora_family.reset_index()
unique_agora_family= unique_agora_family.rename(columns={0: "family_sort"})       #rename column
agora_family_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_agora_family['family_sort']:  
    family_name_key = str(row)
    temp_df = agora[agora["family"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    agora_family_models_dict[family_name_key] = templist
del unique_agora_family, family_name_key, row, templist
#agora genus    
unique_agora_genus = agora.genus.unique()
unique_agora_genus = pd.DataFrame(unique_agora_genus)    #import to pandas
unique_agora_genus['model_list'] = np.empty((len(unique_agora_genus), 0)).tolist()  #add an empty 'model list' column
unique_agora_genus.reset_index()
unique_agora_genus= unique_agora_genus.rename(columns={0: "genus_sort"})       #rename column
agora_genus_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_agora_genus['genus_sort']:  
    genus_name_key = str(row)
    temp_df = agora[agora["genus"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    agora_genus_models_dict[genus_name_key] = templist
del unique_agora_genus, genus_name_key, row, templist
#agora genus2
unique_agora_genus2 = agora.genus2.unique()
unique_agora_genus2 = pd.DataFrame(unique_agora_genus2)    #import to pandas
unique_agora_genus2['model_list'] = np.empty((len(unique_agora_genus2), 0)).tolist()  #add an empty 'model list' column
unique_agora_genus2.reset_index()
unique_agora_genus2= unique_agora_genus2.rename(columns={0: "genus2_sort"})       #rename column
agora_genus2_models_dict ={}         #this is where I'll store the key-value pairs for genus2- [list of models]
for row in unique_agora_genus2['genus2_sort']:  
    genus2_name_key = str(row)
    temp_df = agora[agora["genus2"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    agora_genus2_models_dict[genus2_name_key] = templist
del unique_agora_genus2, genus2_name_key, row, templist
#agora species
unique_agora_species = agora.species.unique()
unique_agora_species = pd.DataFrame(unique_agora_species)    #import to pandas
unique_agora_species['model_list'] = np.empty((len(unique_agora_species), 0)).tolist()  #add an empty 'model list' column
unique_agora_species.reset_index()
unique_agora_species= unique_agora_species.rename(columns={0: "species_sort"})       #rename column
agora_species_models_dict ={}         #this is where I'll store the key-value pairs for species- [list of models]
for row in unique_agora_species['species_sort']:  
    species_name_key = str(row)
    temp_df = agora[agora["species"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    agora_species_models_dict[species_name_key] = templist
del unique_agora_species, species_name_key, row, templist

###############step 2b: create dictionaries for each taxonomic level MAMBO
#columns in mambo: Family, Genus_new, Genus_extracted, Species_extracted
#mambo family_
unique_mambo_family = mambo.Family.unique()
unique_mambo_family = pd.DataFrame(unique_mambo_family)    #import to pandas
unique_mambo_family['model_list'] = np.empty((len(unique_mambo_family), 0)).tolist()  #add an empty 'model list' column
unique_mambo_family.reset_index()
unique_mambo_family= unique_mambo_family.rename(columns={0: "family_sort"})       #rename column
unique_mambo_family.dropna(inplace=True) #extra for mambo, for some reason
family_specific_mambo = mambo.dropna(subset = ["Family"])
mambo_family_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_mambo_family['family_sort']:  
    family_name_key = str(row)
    temp_df = family_specific_mambo[family_specific_mambo["Family"].str.contains(str(row))]        #temp_df should now contain only mambo rows that match that row in unique_mambo_genera2 
    templist = []
    for row in temp_df['model_id']:
        templist.append(str(row))
    del temp_df
    mambo_family_models_dict[family_name_key] = templist
del unique_mambo_family, family_name_key, row, templist, family_specific_mambo
#mambo genus_new
unique_mambo_Genus_new = mambo.Genus_new.unique()
unique_mambo_Genus_new = pd.DataFrame(unique_mambo_Genus_new)    #import to pandas
unique_mambo_Genus_new['model_list'] = np.empty((len(unique_mambo_Genus_new), 0)).tolist()  #add an empty 'model list' column
unique_mambo_Genus_new.reset_index()
unique_mambo_Genus_new= unique_mambo_Genus_new.rename(columns={0: "Genus_new_sort"})       #rename column
unique_mambo_Genus_new.dropna(inplace=True) #extra for mambo, for some reason
Genus_new_specific_mambo = mambo.dropna(subset = ["Genus_new"])
mambo_Genus_new_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_mambo_Genus_new['Genus_new_sort']:  
    Genus_new_name_key = str(row)
    temp_df = Genus_new_specific_mambo[Genus_new_specific_mambo["Genus_new"].str.contains(str(row))]        #temp_df should now contain only mambo rows that match that row in unique_mambo_genera2 
    templist = []
    for row in temp_df['model_id']:
        templist.append(str(row))
    del temp_df
    mambo_Genus_new_models_dict[Genus_new_name_key] = templist
del unique_mambo_Genus_new, Genus_new_name_key, row, templist, Genus_new_specific_mambo
#mambo genus_extracted
unique_mambo_Genus_extracted = mambo.Genus_extracted.unique()
unique_mambo_Genus_extracted = pd.DataFrame(unique_mambo_Genus_extracted)    #import to pandas
unique_mambo_Genus_extracted['model_list'] = np.empty((len(unique_mambo_Genus_extracted), 0)).tolist()  #add an empty 'model list' column
unique_mambo_Genus_extracted.reset_index()
unique_mambo_Genus_extracted= unique_mambo_Genus_extracted.rename(columns={0: "Genus_extracted_sort"})       #rename column
unique_mambo_Genus_extracted.dropna(inplace=True) #extra for mambo, for some reason
Genus_extracted_specific_mambo = mambo.dropna(subset = ["Genus_extracted"])
mambo_Genus_extracted_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_mambo_Genus_extracted['Genus_extracted_sort']:  
    Genus_extracted_name_key = str(row)
    temp_df = Genus_extracted_specific_mambo[Genus_extracted_specific_mambo["Genus_extracted"].str.contains(str(row))]        #temp_df should now contain only mambo rows that match that row in unique_mambo_genera2 
    templist = []
    for row in temp_df['model_id']:
        templist.append(str(row))
    del temp_df
    mambo_Genus_extracted_models_dict[Genus_extracted_name_key] = templist
del unique_mambo_Genus_extracted, Genus_extracted_name_key, row, templist, Genus_extracted_specific_mambo
#mambo species_extracted
unique_mambo_Species_extracted = mambo.Species_extracted.unique()
unique_mambo_Species_extracted = pd.DataFrame(unique_mambo_Species_extracted)    #import to pandas
unique_mambo_Species_extracted['model_list'] = np.empty((len(unique_mambo_Species_extracted), 0)).tolist()  #add an empty 'model list' column
unique_mambo_Species_extracted.reset_index()
unique_mambo_Species_extracted= unique_mambo_Species_extracted.rename(columns={0: "Species_extracted_sort"})       #rename column
unique_mambo_Species_extracted.dropna(inplace=True) #extra for mambo, for some reason
Species_extracted_specific_mambo = mambo.dropna(subset = ["Species_extracted"])
mambo_Species_extracted_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_mambo_Species_extracted['Species_extracted_sort']:  
    Species_extracted_name_key = str(row)
    temp_df = Species_extracted_specific_mambo[Species_extracted_specific_mambo["Species_extracted"].str.contains(str(row))]        #temp_df should now contain only mambo rows that match that row in unique_mambo_genera2 
    templist = []
    for row in temp_df['model_id']:
        templist.append(str(row))
    del temp_df
    mambo_Species_extracted_models_dict[Species_extracted_name_key] = templist
del unique_mambo_Species_extracted, Species_extracted_name_key, row, templist, Species_extracted_specific_mambo

###############step 3: create new otu tables based on 'id to species level', 'id to genus level', 'id to family level' without replacement
#note that this section is generalized; don't need to repeat for MAMBO.
# "fgs_level_sort", "fg_level_sort", "gs_level_sort", "fgs2_level_sort", "fg2_level_sort", "gs2_level_sort"
otus_plus_tax = copy.deepcopy(otus_clean)
#now bin by tax level
otus_to_species_level = otus_plus_tax[otus_plus_tax['Species'].notnull()]   #hey, look at that
otus_to_genus_level = otus_plus_tax[otus_plus_tax['Genus'].notnull() & otus_plus_tax['Species'].isnull()]   #you may award me the genuis medal now
otus_to_family_level = otus_plus_tax[otus_plus_tax['Family'].notnull() & otus_plus_tax['Genus'].isnull() & otus_plus_tax['Species'].isnull()]   
#all right, now things are gonna get a bit weird
#can I get away with iterating through otus_plus_tax and just... adding whatever models fit? And keeping the col with the biggest # models?

#############step 4a: AGORA: map models and reorganize columns
#First, do AGORA; then do MAMBO
otus_plus_tax_plus_models = copy.deepcopy(otus_plus_tax)
otus_plus_tax_plus_models["agora_family_matched_model_list"] = otus_plus_tax_plus_models["Family"].map(agora_family_models_dict)
otus_plus_tax_plus_models["agora_genus_matched_model_list"] = otus_plus_tax_plus_models["Genus"].map(agora_genus_models_dict)
otus_plus_tax_plus_models["agora_genus2_matched_model_list"] = otus_plus_tax_plus_models["Genus"].map(agora_genus2_models_dict)
otus_plus_tax_plus_models["agora_species_matched_model_list"] = otus_plus_tax_plus_models["Species"].map(agora_species_models_dict)

#first, fill in nan values with blank lists.
#make a list of the columns which need nans replaced with blank lists
list_of_cols = ["agora_family_matched_model_list", "agora_genus_matched_model_list", "agora_genus2_matched_model_list", "agora_species_matched_model_list"]

#for each of these olumns, replace all nans with [] (blank lists)
for col in list_of_cols:
    otus_plus_tax_plus_models.loc[otus_plus_tax_plus_models[col].isnull(),[col]] = \
    otus_plus_tax_plus_models.loc[otus_plus_tax_plus_models[col].isnull(),col].apply(lambda x: [])

#create a new column. all_agora_matched_models, that contains a merged list of all matched models at all levels
otus_plus_tax_plus_models["all_agora_matched_models"] = otus_plus_tax_plus_models["agora_family_matched_model_list"] + \
    otus_plus_tax_plus_models["agora_genus_matched_model_list"] + otus_plus_tax_plus_models["agora_genus2_matched_model_list"] + \
        otus_plus_tax_plus_models["agora_species_matched_model_list"]

#this addition created some duplicates (e.g. where gs and gs2 pulled the same models)    
otus_plus_tax_plus_models['all_agora_matched_models'] = otus_plus_tax_plus_models['all_agora_matched_models'].apply(lambda x: list(set(x)))

#nice. Now, we move onto MAMBO. But first, split this datframe into "has agora model" and "needs mambo model"
#note that we can't use the 'null' thing because there are blank lists in this column
otus_with_agora_models = otus_plus_tax_plus_models[~otus_plus_tax_plus_models['all_agora_matched_models'].str.len().eq(0)]  
otus_without_agora_models = otus_plus_tax_plus_models[otus_plus_tax_plus_models['all_agora_matched_models'].str.len().eq(0)]  

#NEW as of 07.27.21: make sure OTUs aren't duplicated across with_models and without_models
unique_otu_ids = otus_with_agora_models['OTU_ID'].unique().tolist()
#then, filter the without_models dataframe to remove OTU_ids that match this list (i.e. already have models)
otus_without_agora_models= otus_without_agora_models[~otus_without_agora_models.OTU_ID.isin(unique_otu_ids)]

#save AGORA model paired taxa:
otus_with_agora_models.to_csv("07.22.21_Hale2018_OTUs_with_agora_models_single_tax_level_picks.csv")
otus_without_agora_models.to_csv("07.22.21_Hale2018_OTUs_without_agora_models_single_tax_level_picks.csv")


###############step 4b: MAMBO map models and reorganize columns 
#recall that most df management is done above, in step 3.
otus_plus_tax2 = copy.deepcopy(otus_without_agora_models)
otus_plus_tax_plus_models2 = copy.deepcopy(otus_plus_tax2)
otus_plus_tax_plus_models2["mambo_family_matched_model_list"] = otus_plus_tax_plus_models2["Family"].map(mambo_family_models_dict)
otus_plus_tax_plus_models2["mambo_genus_new_matched_model_list"] = otus_plus_tax_plus_models2["Genus"].map(mambo_Genus_new_models_dict)
otus_plus_tax_plus_models2["mambo_genus_extracted_matched_model_list"] = otus_plus_tax_plus_models2["Genus"].map(mambo_Genus_extracted_models_dict)
otus_plus_tax_plus_models2["mambo_species_extracted_matched_model_list"] = otus_plus_tax_plus_models2["Species"].map(mambo_Species_extracted_models_dict)

#make a list of the columns which need nans replaced with blank lists
list_of_cols2 = ["mambo_family_matched_model_list", "mambo_genus_new_matched_model_list", "mambo_genus_extracted_matched_model_list", "mambo_species_extracted_matched_model_list"]

#for each of these olumns, replace all nans with [] (blank lists)
for col in list_of_cols2:
    otus_plus_tax_plus_models2.loc[otus_plus_tax_plus_models2[col].isnull(),[col]] = \
    otus_plus_tax_plus_models2.loc[otus_plus_tax_plus_models2[col].isnull(),col].apply(lambda x: [])

#create a new column. all_agora_matched_models, that contains a merged list of all matched models at all levels
otus_plus_tax_plus_models2["all_mambo_matched_models"] = otus_plus_tax_plus_models2["mambo_family_matched_model_list"] + \
    otus_plus_tax_plus_models2["mambo_genus_new_matched_model_list"] + \
        otus_plus_tax_plus_models2["mambo_genus_extracted_matched_model_list"] + \
            otus_plus_tax_plus_models2["mambo_species_extracted_matched_model_list"] 
            
#this addition created some duplicates (e.g. where gs and gs2 pulled the same models)    
otus_plus_tax_plus_models2['all_mambo_matched_models'] = otus_plus_tax_plus_models2['all_mambo_matched_models'].apply(lambda x: list(set(x)))

#nice. Now, we move onto MAMBO. But first, split this datframe into "has agora model" and "needs mambo model"
#note that we can't use the 'null' thing because there are blank lists in this column
otus_with_mambo_models = otus_plus_tax_plus_models2[~otus_plus_tax_plus_models2['all_mambo_matched_models'].str.len().eq(0)]  
otus_without_mambo_models = otus_plus_tax_plus_models2[otus_plus_tax_plus_models2['all_mambo_matched_models'].str.len().eq(0)]  

#NEW as of 07.27.21: make sure OTUs aren't duplicated across with_models and without_models
unique_otu_ids2 = otus_with_mambo_models['OTU_ID'].unique().tolist()
#then, filter the without_models dataframe to remove OTU_ids that match this list (i.e. already have models)
otus_without_mambo_models= otus_without_mambo_models[~otus_without_mambo_models.OTU_ID.isin(unique_otu_ids2)]

#save MAMBO model paired taxa:
otus_with_mambo_models.to_csv("07.22.21_Hale2018_OTUs_with_mambo_models_single_tax_level_picks.csv")
otus_without_mambo_models.to_csv("07.22.21_Hale2018_OTUs_without_mambo_or_agora_models_single_tax_level_picks.csv")

#########UPDATES: have found a few instances of where OTU family names are similar to, but not quite matched to, models
#how to best fix these...
#use otus_without_mambo_models: add a new column called hand_assigned_family, fill in model-matched families, and redo
#family-level_picking for AGORA and MAMBO.

#first, create a dictionary of old-new family names: keys are old, values are new.
model_matched_families_dict = {'Bacteroidales_Incertae_Sedis': 'Bacillales Incertae Sedis XI',\
                          'Clostridiaceae_1': 'Clostridiaceae', 'Clostridiales_vadinBB60_group': 'Clostridiales',\
                              'Family_XIII': 'Clostridiales Incertae Sedis XIII'}
#then, map this dict to otus_without_models
hand_matched_models = copy.deepcopy(otus_without_mambo_models)    
hand_matched_models["hand_assigned_family"] = hand_matched_models["Family"].map(model_matched_families_dict)
#nice, that's going to clear up a lot of issues.

#now, redo AGORA and MAMBO family matching.
hand_matched_models["agora_family_matched_model_list_hand_assigned"] = hand_matched_models["hand_assigned_family"].map(agora_family_models_dict)
otus_with_hand_matched_agora_models = hand_matched_models.dropna(subset=['agora_family_matched_model_list_hand_assigned'])  
otus_without_hand_matched_agora_models = hand_matched_models[hand_matched_models['agora_family_matched_model_list_hand_assigned'].isnull()]
#(note that this is a bit different than code above for splitting into matched vs unmatched- whatever works.)
#NEW as of 07.27.21: make sure OTUs aren't duplicated across with_models and without_models
unique_otu_ids3 = otus_with_hand_matched_agora_models['OTU_ID'].unique().tolist()
#then, filter the without_models dataframe to remove OTU_ids that match this list (i.e. already have models)
otus_without_hand_matched_agora_models= otus_without_hand_matched_agora_models[~otus_without_hand_matched_agora_models.OTU_ID.isin(unique_otu_ids3)]

otus_with_hand_matched_agora_models.to_csv("07.26.21_Hale2018_OTUs_with_agora_models_single_tax_level_hand_picks.csv")

hand_matched_models2= copy.deepcopy(otus_without_hand_matched_agora_models)
hand_matched_models2["mambo_family_matched_model_list_hand_assigned"] = hand_matched_models["hand_assigned_family"].map(mambo_family_models_dict)
otus_with_hand_matched_mambo_models = hand_matched_models2.dropna(subset=['mambo_family_matched_model_list_hand_assigned'])  
otus_without_hand_matched_mambo_models = hand_matched_models2[hand_matched_models2['mambo_family_matched_model_list_hand_assigned'].isnull()]
#(note that this is a bit different than code above for splitting into matched vs unmatched- whatever works.)
#NEW as of 07.27.21: make sure OTUs aren't duplicated across with_models and without_models
unique_otu_ids4 = otus_with_hand_matched_mambo_models['OTU_ID'].unique().tolist()
#then, filter the without_models dataframe to remove OTU_ids that match this list (i.e. already have models)
otus_without_hand_matched_mambo_models= otus_without_hand_matched_mambo_models[~otus_without_hand_matched_mambo_models.OTU_ID.isin(unique_otu_ids4)]

otus_with_hand_matched_mambo_models.to_csv("07.26.21_Hale2018_OTUs_with_mambo_models_single_tax_level_hand_picks.csv")
otus_without_hand_matched_mambo_models.to_csv("07.26.21_Hale2018_OTUs_without_mambo_models_single_tax_level_hand_picks.csv")













