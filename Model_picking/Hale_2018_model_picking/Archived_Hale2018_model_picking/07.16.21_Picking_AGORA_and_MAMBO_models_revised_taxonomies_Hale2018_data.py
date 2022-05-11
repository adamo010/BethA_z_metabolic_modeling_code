#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 10:08:32 2021

@author: adamo010
"""

#The goal here is to match AGORA and MAMBO models to OTUs from Hale 2018 data.
#Specifically, I want to pick models with parts of taxonomies, rather than the whole thing.
#e.g. pick based on genus-species only; this is to account for changes in taxonomies between
#when the AGORA/MAMBO model taxonomies were generated, and when Hale2018 OTU table was generated.

#all this is based on code from 07.15.21_Picking_MAMBO_models_indiv_specific_comms_Hale2018_data.py
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
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Hale_2018_model_picking")

agora = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/AGORA_model_phylogeny_table.csv")
mambo = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/07.23.20_MAMBO_model_taxonomy_assigned_plus_Treponema.csv")
otus = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Hale_2018_model_picking/07.15.21_Hale2018_OTUs_without_MAMBO_models_V1_indiv_specific_comms.csv") #this is the new contaminant-free table

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

###############step 2: create dictionaries for each taxonomic level AGORA
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

###############step 3: create new otu tables based on 'id to species level', 'id to genus level', 'id to family level' without replacement
# "fgs_level_sort", "fg_level_sort", "gs_level_sort", "fgs2_level_sort", "fg2_level_sort", "gs2_level_sort"
otus_plus_tax = copy.deepcopy(otus)
otus_plus_tax = otus_plus_tax.drop(columns={"Unnamed: 0", "species_matched_model_list", "genus_matched_model_list", "family_matched_model_list"})
#add new sorting column
otus_plus_tax["agora_fgs_level_sort"] = otus_plus_tax['Family'] + "; " + otus_plus_tax['Genus'] + "; " + otus_plus_tax['Species']
otus_plus_tax["agora_fg_level_sort"] = otus_plus_tax['Family'] + "; " + otus_plus_tax['Genus']
otus_plus_tax["agora_gs_level_sort"] = otus_plus_tax['Genus'] + "; " + otus_plus_tax['Species']
#now bin by tax level
otus_to_species_level = otus_plus_tax[otus_plus_tax['Species'].notnull()]   #hey, look at that
otus_to_genus_level = otus_plus_tax[otus_plus_tax['Genus'].notnull() & otus_plus_tax['Species'].isnull()]   #you may award me the genuis medal now
otus_to_family_level = otus_plus_tax[otus_plus_tax['Family'].notnull() & otus_plus_tax['Genus'].isnull() & otus_plus_tax['Species'].isnull()]   
#all right, now things are gonna get a bit weird
#can I get away with iterating through otus_plus_tax and just... adding whatever models fit? And keeping the col with the biggest # models?
otus_plus_tax_plus_models = copy.deepcopy(otus_plus_tax)
otus_plus_tax_plus_models["agora_fgs_matched_model_list"] = otus_plus_tax_plus_models["agora_fgs_level_sort"].map(agora_fgs_models_dict)
otus_plus_tax_plus_models["agora_fgs2_matched_model_list"] = otus_plus_tax_plus_models["agora_fgs_level_sort"].map(agora_fgs2_models_dict)
otus_plus_tax_plus_models["agora_fg_matched_model_list"] = otus_plus_tax_plus_models["agora_fg_level_sort"].map(agora_fg_models_dict)
otus_plus_tax_plus_models["agora_fg2_matched_model_list"] = otus_plus_tax_plus_models["agora_fg_level_sort"].map(agora_fg2_models_dict)
otus_plus_tax_plus_models["agora_gs_matched_model_list"] = otus_plus_tax_plus_models["agora_gs_level_sort"].map(agora_gs_models_dict)
otus_plus_tax_plus_models["agora_gs2_matched_model_list"] = otus_plus_tax_plus_models["agora_gs_level_sort"].map(agora_gs2_models_dict)

#what's next? For each row, select the fgs to gs2 column with the longest list.
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
otus_plus_tax_plus_models.to_csv("07.27.21_tracking_down_missing_OTUs.csv")

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

#save AGORA model paired taxa:
otus_with_agora_models.to_csv("07.19.21_Hale2018_OTUs_with_agora_models_round2_picks_indiv_specific_comms.csv")
otus_without_agora_models.to_csv("07.19.21_Hale2018_OTUs_without_agora_models_round2_picks_indiv_specific_comms.csv")

###############step 4: create dictionaries for each taxonomic level MAMBO
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

###############step 5: create new otu tables based on 'id to species level', 'id to genus level', 'id to family level' without replacement
#note that we're going to use otus_without_agora_models; don't want to re-assign models to OTUs which already have AGORA models
###############step 3: create new otu tables based on 'id to species level', 'id to genus level', 'id to family level' without replacement
# "fgs_level_sort", "fg_level_sort", "gs_level_sort", "fgs2_level_sort", "fg2_level_sort", "gs2_level_sort"
otus_plus_tax2 = copy.deepcopy(otus_without_agora_models)
#otus_plus_tax = otus_plus_tax.drop(columns={"Unnamed: 0", "species_matched_model_list", "genus_matched_model_list", "family_matched_model_list"})
#add new sorting column: can we double this up for AGORA?
otus_plus_tax2["mambo_fgs_level_sort"] = otus_plus_tax2['Family'] + "; " + otus_plus_tax2['Genus'] + "; " + otus_plus_tax2['Species']
otus_plus_tax2["mambo_fg_level_sort"] = otus_plus_tax2['Family'] + "; " + otus_plus_tax['Genus']
otus_plus_tax2["mambo_gs_level_sort"] = otus_plus_tax2['Genus'] + "; " + otus_plus_tax['Species']
#now bin by tax level: already done for AGORA
#otus_to_species_level = otus_plus_tax[otus_plus_tax['Species'].notnull()]   #hey, look at that
#otus_to_genus_level = otus_plus_tax[otus_plus_tax['Genus'].notnull() & otus_plus_tax['Species'].isnull()]   #you may award me the genuis medal now
#otus_to_family_level = otus_plus_tax[otus_plus_tax['Family'].notnull() & otus_plus_tax['Genus'].isnull() & otus_plus_tax['Species'].isnull()]   
#all right, now things are gonna get a bit weird
#can I get away with iterating through otus_plus_tax and just... adding whatever models fit? And keeping the col with the biggest # models?

#all right, now things are gonna get a bit weird
#can I get away with iterating through otus_plus_tax and just... adding whatever models fit? And keeping the col with the biggest # models?
otus_plus_tax_plus_models2 = copy.deepcopy(otus_plus_tax2)
otus_plus_tax_plus_models2["mambo_fgs_matched_model_list"] = otus_plus_tax_plus_models2["mambo_fgs_level_sort"].map(mambo_fgs_models_dict)
otus_plus_tax_plus_models2["mambo_fgs2_matched_model_list"] = otus_plus_tax_plus_models2["mambo_fgs_level_sort"].map(mambo_fgs2_models_dict)
otus_plus_tax_plus_models2["mambo_fg_matched_model_list"] = otus_plus_tax_plus_models2["mambo_fg_level_sort"].map(mambo_fg_models_dict)
otus_plus_tax_plus_models2["mambo_fg2_matched_model_list"] = otus_plus_tax_plus_models2["mambo_fg_level_sort"].map(mambo_fg2_models_dict)
otus_plus_tax_plus_models2["mambo_gs_matched_model_list"] = otus_plus_tax_plus_models2["mambo_gs_level_sort"].map(mambo_gs_models_dict)
otus_plus_tax_plus_models2["mambo_gs2_matched_model_list"] = otus_plus_tax_plus_models2["mambo_gs_level_sort"].map(mambo_gs2_models_dict)

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

#save MAMBO model paired taxa:
otus_with_mambo_models.to_csv("07.19.21_Hale2018_OTUs_with_mambo_models_round2_picks_indiv_specific_comms.csv")
otus_without_mambo_models.to_csv("07.19.21_Hale2018_OTUs_without_mambo_or_agora_models_round2_picks_indiv_specific_comms.csv")




