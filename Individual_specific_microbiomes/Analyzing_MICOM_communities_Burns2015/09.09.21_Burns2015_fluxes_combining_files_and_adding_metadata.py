#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 15:21:01 2021

@author: adamo010
"""
##newerer starting point: Date_TBA_combining_flux_values_from_MSI_runs.py; adapting this from Burns to Hale data
#(and from 09.02.21_Hale2018_fluxes_combining_files_and_adding_metadata.py)
#newer starting point: 06.16.21_combinding_TDpointonefive_flux_values_from_MSI_runs.py
#starting point: 05.26.21_combining_TDpointone_inputoutput_values_from_MSI_runs.py

#goal: take all fluxes from MSI run, and organize them into something that can be correlated.

import numpy as np
import cobra
import csv
import subprocess
import os
import pandas as pd
import fnmatch
import shutil
import glob
import re
import copy

######step 1: import files
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Common_growth_medium_runs/June2021_rerunning_MICOM_on_MSI_with_fluxes_on/Sept2021_V30_output")
#NOTE: change this when actual files come in. 

filelist = glob.glob('*_V30_Burns_fluxes_on_fluxes.csv') #why haven't I been doing this forever? this creates a list of files with the specified string

sample_names_list = []
for file in filelist:
    sample_name = re.sub("_V30_Burns_fluxes_on_fluxes.csv", "", file) #delete extra crap on string.  
    sample_names_list.append(sample_name)
del file, sample_name

df_list = []
for file in filelist:
    df= pd.read_csv(file)
    df_list.append(df)
del file, df

######NEW step 2: transpose and take max flux value across all compartments (taxa)
df_list2 = []
for df in df_list:
    df=df.set_index('compartment') #set row names to flux pathway names
    df=df.drop("medium", axis=0) #drop the medium column (not exactly sure what it means in terms of flux anyway)
    df_transposed = df.transpose() #flipperoo so flux pathway names are row names, compartments (taxa) are column names
    maxValues = df_transposed.max(axis = 1) #create a series that represents the max value for each flux pathway across compartments (taxa)
    df_transposed['max_flux'] = maxValues #append this series as a column called max_flux. 
    df_slimmed = df_transposed[df_transposed.max_flux != 0] #drop all rows (flux pathways) where maxValue is 0 or less (i.e. no flux)
    df_maxonly = df_slimmed[['max_flux']]
    df_list2.append(df_maxonly)
del df, df_transposed, maxValues, df_slimmed, df_maxonly 

df_dict = dict(zip(sample_names_list, df_list2)) #add back later, for transposed DFs

######step 3: get a list of all unique flux pathway values in each file
all_fluxpath_values = []
for df in df_list2:
    df = df.reset_index()
    fluxpath_list = df["index"].tolist()
    all_fluxpath_values.extend(fluxpath_list)
del df, fluxpath_list

unique_fluxpath_values = set(all_fluxpath_values) 

#CLEANUP
del all_fluxpath_values, df_list, filelist, sample_names_list

######step 4: create a dataframe to append other dataframes to.
combo_fluxpath_df = pd.DataFrame(unique_fluxpath_values)
#combo_fluxpath_df.reset_index(inplace=True)
combo_fluxpath_df.rename(columns={0: "abbreviated_fluxpath_names"}, inplace=True)
#note here that we need both model_metabolite_names (contains EX_ and _m) and abbreviated_metabolite_names (ends trimmed)
#the former is for merging in data; the latter is for assigning longform names

######step 5: append all input_output files to combo_df
for key, value in df_dict.items():
    temp_dict= dict(zip(value.index, value.max_flux)) #convert df to dictionary (makes it easier to merge)
    combo_fluxpath_df[str(key)] = combo_fluxpath_df['abbreviated_fluxpath_names'].map(temp_dict) #append dictionary values as column values named by key
del temp_dict, key, value    

######step 6: rename fluxes into something useful.
#from 10.12.20_community_medium_data_processing.py
#step 5a: import metabolites key. NEW KEY
fluxpath_key = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Common_growth_medium_runs/May_June2021_MICOM_flux_and_pathway_analyses/Microbiome_flux_key_full_June2021.csv", sep=',')#there are a lot of extra columns here: ony really want the first 2
fluxpath_key_slimmed = fluxpath_key.loc[:, 'abbreviation':'subsystem']    #use : to select full column
del fluxpath_key

######step 7: transpose and add sample metadata.
#first, transpose and convert to longform
combo_fluxpath_df_flipped = combo_fluxpath_df.set_index("abbreviated_fluxpath_names").T #set index and tranpose.
combo_fluxpath_df_flipped.reset_index(inplace=True)
combo_fluxpath_df_flipped_longform = combo_fluxpath_df_flipped.melt(id_vars=["index"], var_name="fluxpath_name", value_name="fluxpath_amount")
#here, "index" is columns to keep (id_vars), "fluxpath_name" is category to merge under (var_name),
#and "fluxpath_amount" is name of column with category values (value_name)    
combo_fluxpath_df_flipped_longform = combo_fluxpath_df_flipped_longform.rename(columns={"index": "Sample_ID"}) #name sample ID column
#then, add metabolite metadata
combo_fluxpath_plus_metab_metadata = pd.merge(left=combo_fluxpath_df_flipped_longform, right=fluxpath_key_slimmed, how='left', 
                                          left_on='fluxpath_name', right_on='abbreviation')
combo_fluxpath_plus_metab_metadata = combo_fluxpath_plus_metab_metadata.drop(columns=["abbreviation"])
combo_fluxpath_plus_metab_metadata = combo_fluxpath_plus_metab_metadata.rename(columns={"description": "fluxpath_description",
                                                                                        "formula": "fluxpath_formula",
                                                                                        "subsystem": "fluxpath_subsystem"}) #name sample ID column

#then, add sample metadata
metadata = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Burns_2015_sample_metadata.csv")#wtf, why are therea bunch of extra rows here?
metadata.dropna(subset = ["Patient_Blind_ID"], inplace=True)
combo_fluxpath_df_plus_metab_sample_metadata = pd.merge(left=combo_fluxpath_plus_metab_metadata, right=metadata, 
                                                        how='left', left_on='Sample_ID', right_on='SampleID') #merge metadata
combo_fluxpath_df_plus_metab_sample_metadata.drop(['SampleID', 'Site', 'MSI_status', 'Stage'], axis=1, inplace=True)#drop superfluous columns
#export!
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/") #analysis directory
combo_fluxpath_df_plus_metab_sample_metadata.to_csv("09.09.21_Burns2015_fluxes_with_metadata_longform.csv")































