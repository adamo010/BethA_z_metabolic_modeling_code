#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 11:01:31 2021

@author: adamo010
"""
#The goal here is to match AGORA and MAMBO models to OTUs from Muehlbauer data.
#Specifically, I want to pick models with parts of taxonomies, rather than the whole thing.
#e.g. pick based on genus-species only; this is to account for changes in taxonomies between
#when the AGORA/MAMBO model taxonomies were generated, and when Hale2018 OTU table was generated.
#also, pick based on single tax levels after that. 

#NEW in 10.19 version3- a new contaminant table
#NEW in 10.15 version2- a new contaminant table
#10.15 version based on combining 07.20.21_Picking_AGORA_and_MAMBO_models_single_tax_levels_Hale2018_data.py 
#and 07.19.21_Picking_AGORA_and_MAMBO_models_revised_taxonomies_Hale2018_data_clean.py

#older versions is based on code from 07.15.21_Picking_MAMBO_models_indiv_specific_comms_Hale2018_data.py
#and 07.13.21_Picking_AGORA_models_indiv_specific_comms_Hale2018_data.py

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
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Muehlbauer_2021_model_picking/")

agora = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/AGORA_model_phylogeny_table.csv")
mambo = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/07.23.20_MAMBO_model_taxonomy_assigned_plus_Treponema.csv")
otus = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Muehlbauer_2021_model_picking/12.02.21_Muehlbauer2021_OTUs_without_MAMBO_models_V1_indiv_specific_comms.csv") 
#for OTU table, V1 and V2 are identical, so that's fine.

#DISCLAIMER AT THE BEGINNING! Niccolai2020 does not have species level IDs. So anything involving species from previous versions is

###############step 1: add columns to agora and mambo that can be matched to otus
#for the 07.16.21 version, I am setting up different columns. Don't want to match whole phylogeny; just parts
#using the columns created in the contaminant filtering pipeline
adapted_agora = copy.deepcopy(agora)  #copy to keep everything clean
adapted_agora["fgs_level_sort"] = adapted_agora['family'] + "; " + adapted_agora['genus'] + "; " + adapted_agora['species'] 
adapted_agora["fg_level_sort"] = adapted_agora['family'] + "; " + adapted_agora['genus']
adapted_agora["gs_level_sort"] = adapted_agora['genus'] + "; " + adapted_agora['species']
adapted_agora["fgs2_level_sort"] = adapted_agora['family'] + "; " + adapted_agora['genus2'] + "; " + adapted_agora['species'] 
adapted_agora["fg2_level_sort"] = adapted_agora['family'] + "; " + adapted_agora['genus2']
adapted_agora["gs2_level_sort"] = adapted_agora['genus2'] + "; " + adapted_agora['species']

#family_genus_species, AND family_genus, AND genus_species, and each individually. 
adapted_mambo = copy.deepcopy(mambo)  #copy to keep everything clean
adapted_mambo["fgs_level_sort"] = adapted_mambo['Family'] + "; " + adapted_mambo['Genus_extracted'] + "; " + adapted_mambo['Species_extracted'] 
adapted_mambo["fg_level_sort"] = adapted_mambo['Family'] + "; " + adapted_mambo['Genus_extracted']
adapted_mambo["gs_level_sort"] = adapted_mambo['Genus_extracted'] + "; " + adapted_mambo['Species_extracted']
adapted_mambo["fgs2_level_sort"] = adapted_mambo['Family'] + "; " + adapted_mambo['Genus_new'] + "; " + adapted_mambo['Species_extracted'] 
adapted_mambo["fg2_level_sort"] = adapted_mambo['Family'] + "; " + adapted_mambo['Genus_new']
adapted_mambo["gs2_level_sort"] = adapted_mambo['Genus_new'] + "; " + adapted_mambo['Species_extracted']
adapted_mambo = adapted_mambo.replace(np.nan, "", regex=True)

###############step 2a: create dictionaries for each taxonomic level AGORA
#agora family_genus_species
unique_agora_fgs = adapted_agora.fgs_level_sort.unique()
unique_agora_fgs = pd.DataFrame(unique_agora_fgs)    #import to pandas
unique_agora_fgs['model_list'] = np.empty((len(unique_agora_fgs), 0)).tolist()  #add an empty 'model list' column
unique_agora_fgs.reset_index()
unique_agora_fgs= unique_agora_fgs.rename(columns={0: "fgs_level_sort"})       #rename column
agora_fgs_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_agora_fgs['fgs_level_sort']:  
    fgs_name_key = str(row)
    temp_df = adapted_agora[adapted_agora["fgs_level_sort"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    agora_fgs_models_dict[fgs_name_key] = templist
    del templist, fgs_name_key, row    
del unique_agora_fgs

#agora family_genus
unique_agora_fg = adapted_agora.fg_level_sort.unique()
unique_agora_fg = pd.DataFrame(unique_agora_fg)    #import to pandas
unique_agora_fg['model_list'] = np.empty((len(unique_agora_fg), 0)).tolist()  #add an empty 'model list' column
unique_agora_fg.reset_index()
unique_agora_fg= unique_agora_fg.rename(columns={0: "fg_level_sort"})       #rename column
agora_fg_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_agora_fg['fg_level_sort']:  
    fg_name_key = str(row)
    temp_df = adapted_agora[adapted_agora["fg_level_sort"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    agora_fg_models_dict[fg_name_key] = templist
    del templist, fg_name_key, row    
del unique_agora_fg

#agora genus_species
unique_agora_gs = adapted_agora.gs_level_sort.unique()
unique_agora_gs = pd.DataFrame(unique_agora_gs)    #import to pandas
unique_agora_gs['model_list'] = np.empty((len(unique_agora_gs), 0)).tolist()  #add an empty 'model list' column
unique_agora_gs.reset_index()
unique_agora_gs= unique_agora_gs.rename(columns={0: "gs_level_sort"})       #rename column
agora_gs_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_agora_gs['gs_level_sort']:  
    gs_name_key = str(row)
    temp_df = adapted_agora[adapted_agora["gs_level_sort"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    agora_gs_models_dict[gs_name_key] = templist
    del templist, gs_name_key, row    
del unique_agora_gs

#agora family_genus_species2
unique_agora_fgs2 = adapted_agora.fgs2_level_sort.unique()
unique_agora_fgs2 = pd.DataFrame(unique_agora_fgs2)    #import to pandas
unique_agora_fgs2['model_list'] = np.empty((len(unique_agora_fgs2), 0)).tolist()  #add an empty 'model list' column
unique_agora_fgs2.reset_index()
unique_agora_fgs2= unique_agora_fgs2.rename(columns={0: "fgs2_level_sort"})       #rename column
agora_fgs2_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_agora_fgs2['fgs2_level_sort']:  
    fgs2_name_key = str(row)
    temp_df = adapted_agora[adapted_agora["fgs2_level_sort"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    agora_fgs2_models_dict[fgs2_name_key] = templist
    del templist, fgs2_name_key, row    
del unique_agora_fgs2

#agora family_genus2
unique_agora_fg2 = adapted_agora.fg2_level_sort.unique()
unique_agora_fg2 = pd.DataFrame(unique_agora_fg2)    #import to pandas
unique_agora_fg2['model_list'] = np.empty((len(unique_agora_fg2), 0)).tolist()  #add an empty 'model list' column
unique_agora_fg2.reset_index()
unique_agora_fg2= unique_agora_fg2.rename(columns={0: "fg2_level_sort"})       #rename column
agora_fg2_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_agora_fg2['fg2_level_sort']:  
    fg2_name_key = str(row)
    temp_df = adapted_agora[adapted_agora["fg2_level_sort"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    agora_fg2_models_dict[fg2_name_key] = templist
    del templist, fg2_name_key, row    
del unique_agora_fg2

#agora genus_species2
unique_agora_gs2 = adapted_agora.gs2_level_sort.unique()
unique_agora_gs2 = pd.DataFrame(unique_agora_gs2)    #import to pandas
unique_agora_gs2['model_list'] = np.empty((len(unique_agora_gs2), 0)).tolist()  #add an empty 'model list' column
unique_agora_gs2.reset_index()
unique_agora_gs2= unique_agora_gs2.rename(columns={0: "gs2_level_sort"})       #rename column
agora_gs2_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_agora_gs2['gs2_level_sort']:  
    gs2_name_key = str(row)
    temp_df = adapted_agora[adapted_agora["gs2_level_sort"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    agora_gs2_models_dict[gs2_name_key] = templist
    del templist, gs2_name_key, row    
del unique_agora_gs2

###############step 2b: create dictionaries for each taxonomic level MAMBO
#mambo family_genus_species
unique_mambo_fgs = adapted_mambo.fgs_level_sort.unique()
unique_mambo_fgs = pd.DataFrame(unique_mambo_fgs)    #import to pandas
unique_mambo_fgs['model_list'] = np.empty((len(unique_mambo_fgs), 0)).tolist()  #add an empty 'model list' column
unique_mambo_fgs.reset_index()
unique_mambo_fgs= unique_mambo_fgs.rename(columns={0: "fgs_level_sort"})       #rename column
mambo_fgs_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_mambo_fgs['fgs_level_sort']:  
    fgs_name_key = str(row)
    temp_df = adapted_mambo[adapted_mambo["fgs_level_sort"].str.contains(str(row))]        #temp_df should now contain only mambo rows that match that row in unique_mambo_genera2 
    templist = []
    for row in temp_df['model_file_name']:
        templist.append(str(row))
    del temp_df
    mambo_fgs_models_dict[fgs_name_key] = templist
    del templist, fgs_name_key, row    
del unique_mambo_fgs

#mambo family_genus
unique_mambo_fg = adapted_mambo.fg_level_sort.unique()
unique_mambo_fg = pd.DataFrame(unique_mambo_fg)    #import to pandas
unique_mambo_fg['model_list'] = np.empty((len(unique_mambo_fg), 0)).tolist()  #add an empty 'model list' column
unique_mambo_fg.reset_index()
unique_mambo_fg= unique_mambo_fg.rename(columns={0: "fg_level_sort"})       #rename column
mambo_fg_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_mambo_fg['fg_level_sort']:  
    fg_name_key = str(row)
    temp_df = adapted_mambo[adapted_mambo["fg_level_sort"].str.contains(str(row))]        #temp_df should now contain only mambo rows that match that row in unique_mambo_genera2 
    templist = []
    for row in temp_df['model_file_name']:
        templist.append(str(row))
    del temp_df
    mambo_fg_models_dict[fg_name_key] = templist
    del templist, fg_name_key, row    
del unique_mambo_fg

#mambo genus_species
unique_mambo_gs = adapted_mambo.gs_level_sort.unique()
unique_mambo_gs = pd.DataFrame(unique_mambo_gs)    #import to pandas
unique_mambo_gs['model_list'] = np.empty((len(unique_mambo_gs), 0)).tolist()  #add an empty 'model list' column
unique_mambo_gs.reset_index()
unique_mambo_gs= unique_mambo_gs.rename(columns={0: "gs_level_sort"})       #rename column
mambo_gs_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_mambo_gs['gs_level_sort']:  
    gs_name_key = str(row)
    temp_df = adapted_mambo[adapted_mambo["gs_level_sort"].str.contains(str(row))]        #temp_df should now contain only mambo rows that match that row in unique_mambo_genera2 
    templist = []
    for row in temp_df['model_file_name']:
        templist.append(str(row))
    del temp_df
    mambo_gs_models_dict[gs_name_key] = templist
    del templist, gs_name_key, row    
del unique_mambo_gs

#mambo family_genus_species2
unique_mambo_fgs2 = adapted_mambo.fgs2_level_sort.unique()
unique_mambo_fgs2 = pd.DataFrame(unique_mambo_fgs2)    #import to pandas
unique_mambo_fgs2['model_list'] = np.empty((len(unique_mambo_fgs2), 0)).tolist()  #add an empty 'model list' column
unique_mambo_fgs2.reset_index()
unique_mambo_fgs2= unique_mambo_fgs2.rename(columns={0: "fgs2_level_sort"})       #rename column
mambo_fgs2_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_mambo_fgs2['fgs2_level_sort']:  
    fgs2_name_key = str(row)
    temp_df = adapted_mambo[adapted_mambo["fgs2_level_sort"].str.contains(str(row))]        #temp_df should now contain only mambo rows that match that row in unique_mambo_genera2 
    templist = []
    for row in temp_df['model_file_name']:
        templist.append(str(row))
    del temp_df
    mambo_fgs2_models_dict[fgs2_name_key] = templist
    del templist, fgs2_name_key, row    
del unique_mambo_fgs2

#mambo family_genus2
unique_mambo_fg2 = adapted_mambo.fg2_level_sort.unique()
unique_mambo_fg2 = pd.DataFrame(unique_mambo_fg2)    #import to pandas
unique_mambo_fg2['model_list'] = np.empty((len(unique_mambo_fg2), 0)).tolist()  #add an empty 'model list' column
unique_mambo_fg2.reset_index()
unique_mambo_fg2= unique_mambo_fg2.rename(columns={0: "fg2_level_sort"})       #rename column
mambo_fg2_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_mambo_fg2['fg2_level_sort']:  
    fg2_name_key = str(row)
    temp_df = adapted_mambo[adapted_mambo["fg2_level_sort"].str.contains(str(row))]        #temp_df should now contain only mambo rows that match that row in unique_mambo_genera2 
    templist = []
    for row in temp_df['model_file_name']:
        templist.append(str(row))
    del temp_df
    mambo_fg2_models_dict[fg2_name_key] = templist
    del templist, fg2_name_key, row    
del unique_mambo_fg2

#mambo genus_species2
unique_mambo_gs2 = adapted_mambo.gs2_level_sort.unique()
unique_mambo_gs2 = pd.DataFrame(unique_mambo_gs2)    #import to pandas
unique_mambo_gs2['model_list'] = np.empty((len(unique_mambo_gs2), 0)).tolist()  #add an empty 'model list' column
unique_mambo_gs2.reset_index()
unique_mambo_gs2= unique_mambo_gs2.rename(columns={0: "gs2_level_sort"})       #rename column
mambo_gs2_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_mambo_gs2['gs2_level_sort']:  
    gs2_name_key = str(row)
    temp_df = adapted_mambo[adapted_mambo["gs2_level_sort"].str.contains(str(row))]        #temp_df should now contain only mambo rows that match that row in unique_mambo_genera2 
    templist = []
    for row in temp_df['model_file_name']:
        templist.append(str(row))
    del temp_df
    mambo_gs2_models_dict[gs2_name_key] = templist
    del templist, gs2_name_key, row    
del unique_mambo_gs2


###############step 3: create new otu tables based on 'id to species level', 'id to genus level', 'id to family level' without replacement
#note that this section is generalized; don't need to repeat for MAMBO.
# "fgs_level_sort", "fg_level_sort", "gs_level_sort", "fgs2_level_sort", "fg2_level_sort", "gs2_level_sort"
otus_plus_tax = copy.deepcopy(otus)
otus_plus_tax = otus_plus_tax.drop(columns={"Unnamed: 0", "genus_matched_model_list", "family_matched_model_list"})
#add new sorting column: works for both AGORA and MAMBO
otus_plus_tax["fgs_level_sort"] = otus_plus_tax['family'] + "; " + otus_plus_tax['genus'] + "; " + otus_plus_tax['species']
otus_plus_tax["fg_level_sort"] = otus_plus_tax['family'] + "; " + otus_plus_tax['genus']
otus_plus_tax["gs_level_sort"] = otus_plus_tax['genus'] + "; " + otus_plus_tax['species']
#now bin by tax level
otus_to_species_level = otus_plus_tax[otus_plus_tax['species'].notnull()]   #hey, look at that
otus_to_genus_level = otus_plus_tax[otus_plus_tax['genus'].notnull() & otus_plus_tax['species'].isnull()] #you may award me the genuis medal now
otus_to_family_level = otus_plus_tax[otus_plus_tax['family'].notnull() & otus_plus_tax['genus'].isnull()]
#all right, now things are gonna get a bit weird
#can I get away with iterating through otus_plus_tax and just... adding whatever models fit? And keeping the col with the biggest # models?

#############step 4a: AGORA: map models and reorganize columns
#First, do AGORA; then do MAMBO
otus_plus_tax_plus_models = copy.deepcopy(otus_plus_tax)
otus_plus_tax_plus_models["agora_fgs_matched_model_list"] = otus_plus_tax_plus_models["fgs_level_sort"].map(agora_fgs_models_dict)
otus_plus_tax_plus_models["agora_fgs2_matched_model_list"] = otus_plus_tax_plus_models["fgs_level_sort"].map(agora_fgs2_models_dict)
otus_plus_tax_plus_models["agora_fg_matched_model_list"] = otus_plus_tax_plus_models["fg_level_sort"].map(agora_fg_models_dict)
otus_plus_tax_plus_models["agora_fg2_matched_model_list"] = otus_plus_tax_plus_models["fg_level_sort"].map(agora_fg2_models_dict)
otus_plus_tax_plus_models["agora_gs_matched_model_list"] = otus_plus_tax_plus_models["gs_level_sort"].map(agora_gs_models_dict)
otus_plus_tax_plus_models["agora_gs2_matched_model_list"] = otus_plus_tax_plus_models["gs_level_sort"].map(agora_gs2_models_dict)

#unify all columns with OTUs matched. Do this for AGORA first, then MAMBO (since the output from AGORA is used for MAMBO)
#first, create new column (all_agora_matched_models) that is basically the merger of all the model mapped columns we just created.
#have to do this weird .fillna('') because if columns contain nan (and many of them do, w/o models matched), it wipes out the whole
#merge. Super frustrating.

#first, fill in nan values with blank lists.
#make a list of the columns which need nans replaced with blank lists
list_of_cols = ["agora_fgs2_matched_model_list", "agora_fgs_matched_model_list", "agora_fg_matched_model_list", "agora_gs_matched_model_list", \
                "agora_fgs2_matched_model_list", "agora_fg2_matched_model_list", "agora_gs2_matched_model_list"]

#for each of these olumns, replace all nans with [] (blank lists)
for col in list_of_cols:
    otus_plus_tax_plus_models.loc[otus_plus_tax_plus_models[col].isnull(),[col]] = \
    otus_plus_tax_plus_models.loc[otus_plus_tax_plus_models[col].isnull(),col].apply(lambda x: [])

#create a new column. all_agora_matched_models, that contains a merged list of all matched models at all levels
otus_plus_tax_plus_models["all_agora_matched_models"] = otus_plus_tax_plus_models["agora_gs2_matched_model_list"] + \
    otus_plus_tax_plus_models["agora_gs_matched_model_list"] + otus_plus_tax_plus_models["agora_fg2_matched_model_list"] + \
        otus_plus_tax_plus_models["agora_fg_matched_model_list"] + otus_plus_tax_plus_models["agora_fgs2_matched_model_list"] + \
            otus_plus_tax_plus_models["agora_fgs_matched_model_list"]

#this addition created some duplicates (e.g. where gs and gs2 pulled the same models)    
otus_plus_tax_plus_models['all_agora_matched_models'] = otus_plus_tax_plus_models['all_agora_matched_models'].apply(lambda x: list(set(x)))

#nice. Now, we move onto MAMBO. But first, split this datframe into "has agora model" and "needs mambo model"
#note that we can't use the 'null' thing because there are blank lists in this column
otus_with_agora_models = otus_plus_tax_plus_models[~otus_plus_tax_plus_models['all_agora_matched_models'].str.len().eq(0)]  
otus_without_agora_models = otus_plus_tax_plus_models[otus_plus_tax_plus_models['all_agora_matched_models'].str.len().eq(0)] 

#NEW addition as of 07.27.21: need to remove OTUs from otus_without_AGORA_models that were actually matched, but have the
#wrong taxonomy to actually have a model assigned.
#basically, if OTU_id is in otus_with_agora_models, remove from otus_without_agora_models.
#first, create a list of unique otu_ids from otus_with_agora_models
unique_otu_ids = otus_with_agora_models['OTU_ID'].unique().tolist()
#then, filter the without_models dataframe to remove OTU_ids that match this list (i.e. already have models)
otus_without_agora_models= otus_without_agora_models[~otus_without_agora_models.OTU_ID.isin(unique_otu_ids)]

###############step 4b: MAMBO map models and reorganize columns 
#recall that most df management is done above, in step 3.
otus_plus_tax2 = copy.deepcopy(otus_without_agora_models)
otus_plus_tax_plus_models2 = copy.deepcopy(otus_plus_tax2)
otus_plus_tax_plus_models2["mambo_fgs_matched_model_list"] = otus_plus_tax_plus_models2["fgs_level_sort"].map(mambo_fgs_models_dict)
otus_plus_tax_plus_models2["mambo_fgs2_matched_model_list"] = otus_plus_tax_plus_models2["fgs_level_sort"].map(mambo_fgs2_models_dict)
otus_plus_tax_plus_models2["mambo_fg_matched_model_list"] = otus_plus_tax_plus_models2["fg_level_sort"].map(mambo_fg_models_dict)
otus_plus_tax_plus_models2["mambo_fg2_matched_model_list"] = otus_plus_tax_plus_models2["fg_level_sort"].map(mambo_fg2_models_dict)
otus_plus_tax_plus_models2["mambo_gs_matched_model_list"] = otus_plus_tax_plus_models2["gs_level_sort"].map(mambo_gs_models_dict)
otus_plus_tax_plus_models2["mambo_gs2_matched_model_list"] = otus_plus_tax_plus_models2["gs_level_sort"].map(mambo_gs2_models_dict)

#what's next? For each row, select the fgs to gs2 column with the longest list.
#first, create new column (all_agora_matched_models) that is basically the merger of all the model mapped columns we just created.
#have to do this weird .fillna('') because if columns contain nan (and many of them do, w/o models matched), it wipes out the whole
#merge. Super frustrating.

#first, fill in nan values with blank lists.
#make a list of the columns which need nans replaced with blank lists
list_of_cols2 = ["mambo_fgs2_matched_model_list", "mambo_fgs_matched_model_list", "mambo_fg_matched_model_list", "mambo_gs_matched_model_list", \
                "mambo_fgs2_matched_model_list", "mambo_fg2_matched_model_list", "mambo_gs2_matched_model_list"]

#for each of these olumns, replace all nans with [] (blank lists)
for col in list_of_cols2:
    otus_plus_tax_plus_models2.loc[otus_plus_tax_plus_models2[col].isnull(),[col]] = \
    otus_plus_tax_plus_models2.loc[otus_plus_tax_plus_models2[col].isnull(),col].apply(lambda x: [])

#create a new column. all_agora_matched_models, that contains a merged list of all matched models at all levels
otus_plus_tax_plus_models2["all_mambo_matched_models"] = otus_plus_tax_plus_models2["mambo_gs2_matched_model_list"] + \
    otus_plus_tax_plus_models2["mambo_gs_matched_model_list"] + otus_plus_tax_plus_models2["mambo_fg2_matched_model_list"] + \
        otus_plus_tax_plus_models2["mambo_fg_matched_model_list"] + otus_plus_tax_plus_models2["mambo_fgs2_matched_model_list"] + \
            otus_plus_tax_plus_models2["mambo_fgs_matched_model_list"]

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


###########################Before moving on, let's clear out some jetsam.
del otus, col, list_of_cols, list_of_cols2, agora_fg2_models_dict, agora_fg_models_dict, mambo_fg2_models_dict, mambo_fg_models_dict,\
    otus_plus_tax, otus_plus_tax2, otus_plus_tax_plus_models, otus_plus_tax_plus_models2, unique_otu_ids,\
        unique_otu_ids2, otus_to_genus_level, otus_to_species_level
#probably still some extra stuff there, but it helps.

#NOW, start pulling from 07.20, single-tax level matching
    #again, we are NOT doing species level here.

otus= copy.deepcopy(otus_without_mambo_models) #remake our otus dataframe here
#really need to clear out excess columns- there are too many
sorted(otus)
otus_clean = otus.drop(columns=['agora_fg2_matched_model_list', 'agora_fg_matched_model_list', 'agora_fgs_matched_model_list',\
                                'agora_fgs2_matched_model_list', 'agora_gs_matched_model_list', 'agora_gs2_matched_model_list',\
                                'all_agora_matched_models', 'all_mambo_matched_models', 'fgs_level_sort', 'gs_level_sort', 'fg_level_sort',\
                                    'mambo_fg2_matched_model_list', 'mambo_fg_matched_model_list','mambo_fgs_matched_model_list',\
                                        'mambo_fgs2_matched_model_list', 'mambo_gs_matched_model_list','mambo_gs2_matched_model_list'])

###############step 2a: create dictionaries for each taxonomic level AGORA
 #columns in agora: family, genus, genus2
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
#columns in mambo: Family, Genus_new, Genus_extracted
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
otus_to_species_level = otus_plus_tax[otus_plus_tax['species'].notnull()]   #hey, look at that
otus_to_genus_level = otus_plus_tax[otus_plus_tax['genus'].notnull() & otus_plus_tax['species'].isnull()]   #you may award me the genuis medal now
otus_to_family_level = otus_plus_tax[otus_plus_tax['family'].notnull() & otus_plus_tax['genus'].isnull()]   
#all right, now things are gonna get a bit weird
#can I get away with iterating through otus_plus_tax and just... adding whatever models fit? And keeping the col with the biggest # models?

#############step 4a: AGORA: map models and reorganize columns
#First, do AGORA; then do MAMBO
otus_plus_tax_plus_models = copy.deepcopy(otus_plus_tax)
otus_plus_tax_plus_models["agora_family_matched_model_list"] = otus_plus_tax_plus_models["family"].map(agora_family_models_dict)
otus_plus_tax_plus_models["agora_genus_matched_model_list"] = otus_plus_tax_plus_models["genus"].map(agora_genus_models_dict)
otus_plus_tax_plus_models["agora_genus2_matched_model_list"] = otus_plus_tax_plus_models["genus"].map(agora_genus2_models_dict)
otus_plus_tax_plus_models["agora_species_matched_model_list"] = otus_plus_tax_plus_models["species"].map(agora_species_models_dict)

#first, fill in nan values with blank lists.
#make a list of the columns which need nans replaced with blank lists
list_of_cols = ["agora_family_matched_model_list", "agora_genus_matched_model_list", "agora_genus2_matched_model_list", "agora_species_matched_model_list"]

#for each of these olumns, replace all nans with [] (blank lists)
for col in list_of_cols:
    otus_plus_tax_plus_models.loc[otus_plus_tax_plus_models[col].isnull(),[col]] = \
    otus_plus_tax_plus_models.loc[otus_plus_tax_plus_models[col].isnull(),col].apply(lambda x: [])
del col, list_of_cols

#create a new column. all_agora_matched_models, that contains a merged list of all matched models at all levels
otus_plus_tax_plus_models["all_agora_matched_models"] = otus_plus_tax_plus_models["agora_family_matched_model_list"] + \
    otus_plus_tax_plus_models["agora_genus_matched_model_list"] + otus_plus_tax_plus_models["agora_genus2_matched_model_list"] + \
        otus_plus_tax_plus_models["agora_species_matched_model_list"]

#this addition created some duplicates (e.g. where gs and gs2 pulled the same models)    
otus_plus_tax_plus_models['all_agora_matched_models'] = otus_plus_tax_plus_models['all_agora_matched_models'].apply(lambda x: list(set(x)))

#nice. Now, we move onto MAMBO. But first, split this datframe into "has agora model" and "needs mambo model"
#note that we can't use the 'null' thing because there are blank lists in this column
otus_with_agora_models2 = otus_plus_tax_plus_models[~otus_plus_tax_plus_models['all_agora_matched_models'].str.len().eq(0)]  
otus_without_agora_models2 = otus_plus_tax_plus_models[otus_plus_tax_plus_models['all_agora_matched_models'].str.len().eq(0)]  

#NEW as of 07.27.21: make sure OTUs aren't duplicated across with_models and without_models
unique_otu_ids = otus_with_agora_models['OTU_ID'].unique().tolist()
#then, filter the without_models dataframe to remove OTU_ids that match this list (i.e. already have models)
otus_without_agora_models2= otus_without_agora_models2[~otus_without_agora_models2.OTU_ID.isin(unique_otu_ids)]

#save AGORA model paired taxa:
#otus_with_agora_models.to_csv("07.22.21_Hale2018_OTUs_with_agora_models_single_tax_level_picks.csv")
#otus_without_agora_models.to_csv("07.22.21_Hale2018_OTUs_without_agora_models_single_tax_level_picks.csv")

###############step 4b: MAMBO map models and reorganize columns 
#recall that most df management is done above, in step 3.
otus_plus_tax2 = copy.deepcopy(otus_without_agora_models2)
otus_plus_tax_plus_models2 = copy.deepcopy(otus_plus_tax2)
otus_plus_tax_plus_models2["mambo_family_matched_model_list"] = otus_plus_tax_plus_models2["family"].map(mambo_family_models_dict)
otus_plus_tax_plus_models2["mambo_genus_new_matched_model_list"] = otus_plus_tax_plus_models2["genus"].map(mambo_Genus_new_models_dict)
otus_plus_tax_plus_models2["mambo_genus_extracted_matched_model_list"] = otus_plus_tax_plus_models2["genus"].map(mambo_Genus_extracted_models_dict)
otus_plus_tax_plus_models2["mambo_species_extracted_matched_model_list"] = otus_plus_tax_plus_models2["species"].map(mambo_Species_extracted_models_dict)

#make a list of the columns which need nans replaced with blank lists
list_of_cols2 = ["mambo_family_matched_model_list", "mambo_genus_new_matched_model_list", "mambo_genus_extracted_matched_model_list",\
                 "mambo_species_extracted_matched_model_list"]

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
otus_with_mambo_models2 = otus_plus_tax_plus_models2[~otus_plus_tax_plus_models2['all_mambo_matched_models'].str.len().eq(0)]  
otus_without_mambo_models2 = otus_plus_tax_plus_models2[otus_plus_tax_plus_models2['all_mambo_matched_models'].str.len().eq(0)]  

#NEW as of 07.27.21: make sure OTUs aren't duplicated across with_models and without_models
unique_otu_ids2 = otus_with_mambo_models2['OTU_ID'].unique().tolist()
#then, filter the without_models dataframe to remove OTU_ids that match this list (i.e. already have models)
otus_without_mambo_models2= otus_without_mambo_models2[~otus_without_mambo_models2.OTU_ID.isin(unique_otu_ids2)]

#save our output:
otus_with_agora_models.to_csv("12.02.21_Muehlbauer2021_OTUs_with_agora_models_fg_level_picks.csv")
otus_with_agora_models2.to_csv("12.02.21_Muehlbauer2021_OTUs_with_agora_models_single_tax_level_picks.csv")
otus_without_mambo_models2.to_csv("12.02.21_Muehlbauer2021_OTUs_without_models_after_fg_and_single_tax_level_picking.csv")
#otus_with_mambo_models.to_csv("07.22.21_Hale2018_OTUs_with_mambo_models_single_tax_level_picks.csv")
#otus_without_mambo_models.to_csv("07.22.21_Hale2018_OTUs_without_mambo_or_agora_models_single_tax_level_picks.csv")

#AND NOW FOR SOME HAND MATCHING
#########UPDATES: have found a few instances of where OTU family names are similar to, but not quite matched to, models
#how to best fix these...
#use otus_without_mambo_models: add a new column called hand_assigned_family, fill in model-matched families, and redo
#family-level_picking for AGORA and MAMBO.

#first, create a dictionary of old-new family names: keys are old, values are new.
model_matched_families_dict = {'Clostridiaceae_1': 'Clostridiaceae'}
#then, map this dict to otus_without_models
hand_matched_models = copy.deepcopy(otus_without_mambo_models2)    
hand_matched_models["hand_assigned_family"] = hand_matched_models["family"].map(model_matched_families_dict)
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

otus_with_hand_matched_agora_models.to_csv("12.02.21_Muehlbauer2021_with_agora_models_single_tax_level_hand_picks.csv")

hand_matched_models2= copy.deepcopy(otus_without_hand_matched_agora_models)
hand_matched_models2["mambo_family_matched_model_list_hand_assigned"] = hand_matched_models["hand_assigned_family"].map(mambo_family_models_dict)
otus_with_hand_matched_mambo_models = hand_matched_models2.dropna(subset=['mambo_family_matched_model_list_hand_assigned'])  
otus_without_hand_matched_mambo_models = hand_matched_models2[hand_matched_models2['mambo_family_matched_model_list_hand_assigned'].isnull()]
#(note that this is a bit different than code above for splitting into matched vs unmatched- whatever works.)
#NEW as of 07.27.21: make sure OTUs aren't duplicated across with_models and without_models
unique_otu_ids4 = otus_with_hand_matched_mambo_models['OTU_ID'].unique().tolist()
#then, filter the without_models dataframe to remove OTU_ids that match this list (i.e. already have models)
otus_without_hand_matched_mambo_models= otus_without_hand_matched_mambo_models[~otus_without_hand_matched_mambo_models.OTU_ID.isin(unique_otu_ids4)]

otus_with_hand_matched_mambo_models.to_csv("12.02.21_Muehlbauer2021_OTUs_with_mambo_models_single_tax_level_hand_picks.csv")
otus_without_hand_matched_mambo_models.to_csv("12.02.21_Muehlbauer2021_OTUs_without_mambo_models_single_tax_level_hand_picks.csv")



    
    