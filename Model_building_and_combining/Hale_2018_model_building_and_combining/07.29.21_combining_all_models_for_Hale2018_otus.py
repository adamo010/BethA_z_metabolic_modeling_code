#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 14:46:14 2021

@author: adamo010
"""

#based on 05.26.20_combining_CarveMe_models_for_Burns_otus.py

#the goal of this code is to a) combine all dataframes of agora, mambo, and carveme models,
#b) combine all the models together where multiple models have been assigned to a single OTU
#c) put all models for this particular dataset in a central repository for running later.

import micom
from micom.util import (
    load_model,
    join_models,
    add_var_from_expression,
    adjust_solver_config,
    clean_ids
)
import gurobipy
import scipy
import pickle
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
import ast
from ast import literal_eval

#####Step 1: All AGORA models: import and fix up
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Hale_2018_model_picking/")

#the goal is to get everything down to 102 columns: this is taxonomy_orig_x, OTU_ID, a bunch of sample_ids, then cols 
#kingdom, Phylum, Class, Order, Family, Genus, Species, and all_agora_matched_models

####Fixing up V1
agora_V1 = pd.read_csv("07.13.21_Hale2018_OTUs_with_models_V1_indiv_specific_comms.csv")   
#make a list of the columns which need nans replaced with blank lists
list_of_cols = ["family_matched_model_list", "genus_matched_model_list", "species_matched_model_list"]
#for each of these olumns, replace all nans with [] (blank lists)
for col in list_of_cols:
    agora_V1.loc[agora_V1[col].isnull(),[col]] = \
    agora_V1.loc[agora_V1[col].isnull(),col].apply(lambda x: [])
#now, convert all strings back to lists (literal_eval does this)    
agora_V1['family_matched_model_list'] = agora_V1.family_matched_model_list.apply(lambda x: literal_eval(str(x)))
agora_V1['genus_matched_model_list'] = agora_V1.genus_matched_model_list.apply(lambda x: literal_eval(str(x)))
agora_V1['species_matched_model_list'] = agora_V1.species_matched_model_list.apply(lambda x: literal_eval(str(x)))
#create a merged column with all AGORA models   
agora_V1["all_agora_matched_models"] = agora_V1["family_matched_model_list"] + agora_V1["genus_matched_model_list"] +\
    agora_V1["species_matched_model_list"]

#drop unnecessary columns.
sorted(agora_V1) #prints column names
agora_V1 = agora_V1.drop(columns=["Unnamed: 0",'class_level_sort','clean_species_level_sort','family_level_sort','family_matched_model_list',\
                                 'genus_level_sort','genus_matched_model_list','keep_or_toss','order_level_sort','phylum_level_sort',\
                                     'species_level_sort','taxonomy_orig_y','species_matched_model_list'])  

####Fixing up V2
agora_V2 = pd.read_csv("07.19.21_Hale2018_OTUs_with_agora_models_round2_picks_indiv_specific_comms.csv")
sorted(agora_V2) #prints column names
#for V2, already have an all_models column. so, can delete other matched_model_list columns. 
agora_V2 = agora_V2.drop(columns=["Unnamed: 0", 'keep_or_toss', 'class_level_sort', 'clean_species_level_sort','family_level_sort',\
                                 'genus_level_sort', 'order_level_sort', 'phylum_level_sort', 'species_level_sort',\
                                     "agora_fg2_matched_model_list", "agora_fg_matched_model_list", "agora_fgs2_matched_model_list",\
                                         'agora_fgs_matched_model_list', 'agora_gs2_matched_model_list', 'agora_gs_matched_model_list',\
                                             'fgs_level_sort', 'fg_level_sort', 'gs_level_sort'])
#however, do still need to convert all_agora_models column to list
agora_V2['all_agora_matched_models'] = agora_V2.all_agora_matched_models.apply(lambda x: literal_eval(str(x)))

####Fixing up V3
agora_V3 = pd.read_csv("07.22.21_Hale2018_OTUs_with_agora_models_single_tax_level_picks.csv")
sorted(agora_V3)
#for V3, already have an all_models column. so, can delete other matched_model_list columns. 
agora_V3 = agora_V3.drop(columns=["Unnamed: 0", 'class_level_sort','genus_level_sort', 'order_level_sort', 'phylum_level_sort',\
                                  'species_level_sort', "agora_family_matched_model_list", "agora_genus2_matched_model_list",\
                                      "agora_genus_matched_model_list", 'agora_species_matched_model_list'])
#however, do still need to convert all_agora_models column to list
agora_V3['all_agora_matched_models'] = agora_V3.all_agora_matched_models.apply(lambda x: literal_eval(str(x)))

####Fixing up V4
agora_V4 = pd.read_csv("07.26.21_Hale2018_OTUs_with_agora_models_single_tax_level_hand_picks.csv")
sorted(agora_V4)
#for V3, already have an all_models column. so, can delete other matched_model_list columns. 
agora_V4 = agora_V4.drop(columns=["Unnamed: 0", 'agora_family_matched_model_list','agora_genus2_matched_model_list',\
                                  'agora_genus_matched_model_list', 'agora_species_matched_model_list',\
                                  'all_agora_matched_models', "all_mambo_matched_models", "class_level_sort",\
                                      'genus_level_sort', 'hand_assigned_family','mambo_family_matched_model_list',\
                                          'mambo_genus_extracted_matched_model_list','mambo_genus_new_matched_model_list',\
                                              'mambo_species_extracted_matched_model_list', 'order_level_sort',\
                                                  'phylum_level_sort', 'species_level_sort'])
#note here that the column that should be named all_agora_matched_models is here called agora_family_matched_model_list_hand_assigned.
#need to rename it.
agora_V4.rename(columns = {'agora_family_matched_model_list_hand_assigned':'all_agora_matched_models'}, inplace = True)    
#and, do still need to convert all_agora_models column to list
agora_V4['all_agora_matched_models'] = agora_V4.all_agora_matched_models.apply(lambda x: literal_eval(str(x)))

#now, create a whole merged dataframe of AGORA_matched_models
#all these dataframes should have the same headers, so should be easy
all_agora_merged = pd.concat([agora_V1, agora_V2, agora_V3, agora_V4], axis=0)
#nice. have 2300 rows even. 
#now for the trickier part- merging columns with shared value in OTU_ID
#stackoverflow recommends collapsing (groupby) first, which causes a loss of most columns. Then, re-add back other columns from orig dataframe
all_agora_merged_slimmed = all_agora_merged.groupby('OTU_ID')['all_agora_matched_models'].sum().reset_index()
#this addition created some duplicates (e.g. where gs and gs2 pulled the same models)    
all_agora_merged_slimmed['all_agora_matched_models'] = all_agora_merged_slimmed['all_agora_matched_models'].apply(lambda x: list(set(x)))
#all_agora_merged.to_csv("testoJuly29.csv")
#all_agora_merged_slimmed.to_csv("testoJuly29-2.csv")
#great, seems like that worked. Now, merge back in other data.
#merged_left = pd.merge(left=survey_sub, right=species_sub, how='left', left_on='species_id', right_on='species_id')
all_agora_merged_condensed = pd.merge(left= all_agora_merged_slimmed, right= all_agora_merged, how= 'left',\
                                      left_on="OTU_ID", right_on = "OTU_ID")
#all right, that looks weird. have the same number of rows as all_agora_merged, not all_agora_merged_slimmed.
#an important note here- have to drop_duplicates based on a subset because this command doesn't work on lists    
    #which is annoying, but shouldn't matter.
all_agora_merged_condensed.drop_duplicates(subset="OTU_ID", keep='first', inplace=True)  
all_agora_merged_condensed.drop(columns=["all_agora_matched_models_y"], inplace=True)

#clean up
del agora_V1, agora_V2, agora_V3, agora_V4, all_agora_merged, all_agora_merged_slimmed, col, list_of_cols

#####Step 2: All MAMBO models
####Fixing up V1
mambo_V1 = pd.read_csv('07.15.21_Hale2018_OTUs_with_MAMBO_models_V1_indiv_specific_comms.csv')
#make a list of the columns which need nans replaced with blank lists
list_of_cols = ["family_matched_model_list", "genus_matched_model_list", "species_matched_model_list"]
#for each of these olumns, replace all nans with [] (blank lists)
for col in list_of_cols:
    mambo_V1.loc[mambo_V1[col].isnull(),[col]] = \
    mambo_V1.loc[mambo_V1[col].isnull(),col].apply(lambda x: [])
#now, convert all strings back to lists (literal_eval does this)    
mambo_V1['family_matched_model_list'] = mambo_V1.family_matched_model_list.apply(lambda x: literal_eval(str(x)))
mambo_V1['genus_matched_model_list'] = mambo_V1.genus_matched_model_list.apply(lambda x: literal_eval(str(x)))
mambo_V1['species_matched_model_list'] = mambo_V1.species_matched_model_list.apply(lambda x: literal_eval(str(x)))
#create a merged column with all mambo models   
mambo_V1["all_mambo_matched_models"] = mambo_V1["family_matched_model_list"] + mambo_V1["genus_matched_model_list"] +\
    mambo_V1["species_matched_model_list"]
#drop unnecessary columns.
sorted(mambo_V1) #prints column names
mambo_V1 = mambo_V1.drop(columns=["Unnamed: 0",'class_level_sort','clean_species_level_sort','family_level_sort','family_matched_model_list',\
                                 'genus_level_sort','genus_matched_model_list','keep_or_toss','order_level_sort','phylum_level_sort',\
                                     'species_level_sort','species_matched_model_list'])  
####Fixing up V2
mambo_V2 = pd.read_csv('07.19.21_Hale2018_OTUs_with_mambo_models_round2_picks_indiv_specific_comms.csv')
sorted(mambo_V2) #prints column names
#for V2, already have an all_models column. so, can delete other matched_model_list columns. 
mambo_V2 = mambo_V2.drop(columns=["Unnamed: 0", 'keep_or_toss', 'class_level_sort', 'clean_species_level_sort','family_level_sort',\
                                 'genus_level_sort', 'order_level_sort', 'phylum_level_sort', 'species_level_sort',\
                                     "mambo_fg2_matched_model_list", "mambo_fg_matched_model_list", "mambo_fgs2_matched_model_list",\
                                         'mambo_fgs_matched_model_list', 'mambo_gs2_matched_model_list', 'mambo_gs_matched_model_list',\
                                             "agora_fg2_matched_model_list", "agora_fg_matched_model_list", "agora_fgs2_matched_model_list",\
                                                 'agora_fgs_matched_model_list', 'agora_gs2_matched_model_list', 'agora_gs_matched_model_list',\
                                                     'fgs_level_sort', 'fg_level_sort', 'gs_level_sort', "all_agora_matched_models"])
#however, do still need to convert all_mambo_models column to list
mambo_V2['all_mambo_matched_models'] = mambo_V2.all_mambo_matched_models.apply(lambda x: literal_eval(str(x)))

####Fixing up V3
mambo_V3 = pd.read_csv('07.22.21_Hale2018_OTUs_with_mambo_models_single_tax_level_picks.csv')
sorted(mambo_V3)
#for V3, already have an all_models column. so, can delete other matched_model_list columns. 
mambo_V3 = mambo_V3.drop(columns=["Unnamed: 0", 'agora_family_matched_model_list', 'agora_genus2_matched_model_list','agora_genus_matched_model_list',\
                                  'agora_species_matched_model_list', 'all_agora_matched_models', 'class_level_sort','genus_level_sort',\
                                      "mambo_family_matched_model_list", "mambo_genus_extracted_matched_model_list", "mambo_genus_new_matched_model_list",\
                                          "mambo_species_extracted_matched_model_list", 'order_level_sort', 'phylum_level_sort','species_level_sort']) 
#however, do still need to convert all_mambo_models column to list
mambo_V3['all_mambo_matched_models'] = mambo_V3.all_mambo_matched_models.apply(lambda x: literal_eval(str(x)))

####Fixing up V4
mambo_V4 = pd.read_csv('07.26.21_Hale2018_OTUs_with_mambo_models_single_tax_level_hand_picks.csv')
sorted(mambo_V4)
#for V3, already have an all_models column. so, can delete other matched_model_list columns. 
mambo_V4 = mambo_V4.drop(columns=["Unnamed: 0", 'agora_family_matched_model_list', 'agora_family_matched_model_list_hand_assigned',\
                                  'agora_genus2_matched_model_list', 'agora_genus_matched_model_list', 'agora_species_matched_model_list',\
                                      'all_agora_matched_models', 'all_mambo_matched_models', "class_level_sort",'genus_level_sort',\
                                          'hand_assigned_family','mambo_family_matched_model_list', 'mambo_genus_extracted_matched_model_list',\
                                              'mambo_genus_new_matched_model_list', 'mambo_species_extracted_matched_model_list',\
                                                  'order_level_sort', 'phylum_level_sort','species_level_sort'])
                                              
#note here that the column that should be named all_mambo_matched_models is here called mambo_family_matched_model_list_hand_assigned.
#need to rename it.
mambo_V4.rename(columns = {'mambo_family_matched_model_list_hand_assigned':'all_mambo_matched_models'}, inplace = True)    
#and, do still need to convert all_mambo_models column to list
mambo_V4['all_mambo_matched_models'] = mambo_V4.all_mambo_matched_models.apply(lambda x: literal_eval(str(x)))

#now, create a whole merged dataframe of mambo_matched_models
#all these dataframes should have the same headers, so should be easy
all_mambo_merged = pd.concat([mambo_V1, mambo_V2, mambo_V3, mambo_V4], axis=0)
#nice. have 2300 rows even. 
#now for the trickier part- merging columns with shared value in OTU_ID
#stackoverflow recommends collapsing (groupby) first, which causes a loss of most columns. Then, re-add back other columns from orig dataframe
all_mambo_merged_slimmed = all_mambo_merged.groupby('OTU_ID')['all_mambo_matched_models'].sum().reset_index()
#this addition created some duplicates (e.g. where gs and gs2 pulled the same models)    
all_mambo_merged_slimmed['all_mambo_matched_models'] = all_mambo_merged_slimmed['all_mambo_matched_models'].apply(lambda x: list(set(x)))
#all_mambo_merged.to_csv("testoJuly29.csv")
#all_mambo_merged_slimmed.to_csv("testoJuly29-2.csv")
#great, seems like that worked. Now, merge back in other data.
#merged_left = pd.merge(left=survey_sub, right=species_sub, how='left', left_on='species_id', right_on='species_id')
all_mambo_merged_condensed = pd.merge(left= all_mambo_merged_slimmed, right= all_mambo_merged, how= 'left',\
                                      left_on="OTU_ID", right_on = "OTU_ID")
#all right, that looks weird. have the same number of rows as all_mambo_merged, not all_mambo_merged_slimmed.
#an important note here- have to drop_duplicates based on a subset because this command doesn't work on lists    
    #which is annoying, but shouldn't matter.
all_mambo_merged_condensed.drop_duplicates(subset="OTU_ID", keep='first', inplace=True)  
all_mambo_merged_condensed.drop(columns=["all_mambo_matched_models_y"], inplace=True)

#clean up
del mambo_V1, mambo_V2, mambo_V3, mambo_V4, all_mambo_merged, all_mambo_merged_slimmed, col, list_of_cols

####Step 3: All CarveMe models
carveme = pd.read_csv("07.29.21_Hale2018_OTUs_with_carveme_models.csv")
sorted(carveme)
carveme= carveme.drop(columns=["Unnamed: 0", 'carveme_family_matched_model_list', 'carveme_genus_matched_model_list',\
                               'carveme_species_matched_model_list', 'class_level_sort', 'genus_level_sort', 'hand_assigned_family',\
                                  'order_level_sort', 'phylum_level_sort','species_level_sort'])
carveme['all_carveme_matched_models'] = carveme.all_carveme_matched_models.apply(lambda x: literal_eval(str(x)))

####step 4: import and match OTU_id keys. 
otu_id_key = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Hale2018/Hale_2018_OTUs_with_matched_taxonomy_and_proxy_OTUids.csv")
otu_id_key.drop(columns=["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"], inplace=True)
all_agora_merged_condensed2 = pd.merge(left= all_agora_merged_condensed, right= otu_id_key, how="left", left_on="OTU_ID", right_on="OTU_ID")
all_mambo_merged_condensed2 = pd.merge(left= all_mambo_merged_condensed, right= otu_id_key, how="left", left_on="OTU_ID", right_on="OTU_ID")
carveme2 = pd.merge(left= carveme, right= otu_id_key, how="left", left_on="OTU_ID", right_on="OTU_ID")


####step 5a: import agora model key
agora_model_file_key= pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/agora_model_id_and_file_key.csv",\
                                  index_col=0, squeeze=True).to_dict()

matched_model_id_list = all_agora_merged_condensed2["all_agora_matched_models_x"].to_list() #convert model ID column to list
all_model_files = [] 

#for each model id in the list of model_ids associated with each row (OTU) in the dataframe, match with a model file
#from the agora_model_file_key dictionary
for elem in matched_model_id_list:
    model_file_list = []
    for subelem in elem:
        model_file_string = str(agora_model_file_key[subelem])
        model_file_list.append(model_file_string)
    all_model_files.append(model_file_list)
del elem, subelem, model_file_string    

#now, have a list called all_model_files that can be stuck on as a column in the agora dataframe
all_agora_merged_condensed2['all_agora_model_file_names'] = all_model_files

#cleanup
del all_model_files, agora_model_file_key, matched_model_id_list, model_file_list

#intermediate step: export all data to pickles files, to preserve list format. 
all_agora_merged_condensed2.to_pickle('AGORA_models_that_need_OTU_builds.pickle')
#testo = pd.read_pickle('AGORA_models_that_need_OTU_builds.pickle')
all_mambo_merged_condensed2.to_pickle('MAMBO_models_that_need_OTU_builds.pickle')
carveme2.to_pickle("CarveMe_models_that_need_otu_builds.pickle")


####step 6: merge models. How do we do this?
#####08.04.21: this is taking forever locally; try running on MSI. Everything after this is tesing for MSI run. 


agora_model_file_dict = dict(zip(all_agora_merged_condensed2.OTU_proxy_ID, \
                                 all_agora_merged_condensed2.all_agora_model_file_names.apply(lambda x: literal_eval(str(x)))))    


os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/AGORA_models/")
counter = 0
for key, value in agora_model_file_dict.items():
    filename = str(key)+"_model.xml"
    otu_model = micom.util.join_models(value, str(key))
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/Hale_2018_model_building_and_combining/Hale_2018_otu_specific_models/")
    cobra.io.write_sbml_model(otu_model, filename)
    counter= counter+1
    print(counter)
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/AGORA_models/")







otu_proxy_id_list = all_agora_merged_condensed2["OTU_proxy_ID"].tolist()
combined_model_list = []
for row in all_agora_merged_condensed2["all_agora_model_file_names"]:
    print(type(row))
    otu_model = micom.util.join_models(row, "new_otu_model")
    combined_model_list.append(otu_model)

model_files_dict = dict(zip(otu_proxy_id_list, combined_model_list)) #good thing lists are ordered, that's all I'll say. 

os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/Hale_2018_model_building_and_combining/Hale_2018_otu_specific_models/")

#so, now the trick is to iterate through the dictionary and save each model as an xml file. 
for key, value in model_files_dict.items():
    filename = str(key)+"_model.xml"
    cobra.io.write_sbml_model(value, filename)  #I don't know if this will work, but maybe?



########################################

#rename all items in combined_model_list with the ordered corresponding items in the OTU_proxy_ID column

model_files_dict = dict(zip(otuy_proxy_id_list, combined_model_list)) 





#extract 
agora_model_dict = dict(zip(all_agora_merged_condensed2.OTU_proxy_ID, all_agora_merged_condensed2.all_agora_matched_models_x))
merged_agora_models_dict = {}
for key, value in agora_model_dict.items():
    otu_model = micom.util.join_models(value, str(key))


#hello, I work (from 05.01.20_combining_models_for_Burns_otus.py            
otu_196131_model = micom.util.join_models(["Odoribacter_0.xml", "Odoribacter_1.xml"], "otu_196131_model.xml")
cobra.io.write_sbml_model(otu_196131_model, "otu_196131_model.xml")

#########################Old junk
####step 5: fix the fucking agora models.
#they're listed as model ids rather than .xml file names. pain in the ass.
#iterate over each model_id in each list of each row, and make a new list of lists of file names from model ids. 
#yeah, I know. it fucking sucks
#agora_model_file_names_list=[]
for row in all_agora_merged_condensed2["all_agora_matched_models_x"]:
    otu_model_list = []
    for item in row:
        item2= item.replace(" ", "_")
        item2=item2.replace('-',"_")
        item2=item2.replace('.','')
        item2=item2.replace(',','')
        item2=item2.replace('/','_')
        item2=item2.replace('=_',"")
        item2=item2.replace("__","_")
        item3=item2+".xml"
        item3=item3.replace("_.xml",".xml")
        otu_model_list.append(item3)
    agora_model_file_names_list.append(otu_model_list)    
del item, item2, item3, row
all_agora_merged_condensed2['all_agora_matched_models'] = agora_model_file_names_list
all_agora_merged_condensed2.drop(columns= ["all_agora_matched_models_x"], inplace=True)

#ARGH. Doing this peacemeal sucks. let's do it systemically.
agora_file_names = [v for v in os.listdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/AGORA_models/")] #get a list of all model files
agora_file_names_df = pd.DataFrame(agora_file_names,columns=['model_file_name'])
agora_file_names_df.to_csv("agora_file_names.csv")
#compare this to otu_file_list
model_list_difference = [item for item in agora_file_names if item not in otu_model_list]
