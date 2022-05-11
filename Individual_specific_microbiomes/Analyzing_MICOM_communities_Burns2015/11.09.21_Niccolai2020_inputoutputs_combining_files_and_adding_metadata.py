#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 15:37:51 2021

@author: adamo010
"""
#newerer starting point: Date_TBA_combining_inputoutput_values_from_MSI_runs.py; adapting this from Burns to Hale data
#newer starting point: 05.26.21_combining_TDpointone_inputoutput_values_from_MSI_runs.py
#starting point: V02.04.21_combining_species_specific_GRs_over_TD_values.py AND
#09.02.21_Hale2018_inputoutputs_combining_files_and_adding_metadata.py AND
#09.09.21_Burns2015_inputoutputs_combining_files_and_adding_metadata.py
#tradeoff = 0.2

#goal: take all input/output metabolites from MSI run, and organize them into something that can be correlated.

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

##############################################################################################
###############NEW AND IMPROVED: recent version with improved metabolite key##################
##############################################################################################
#####step 1: import files
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Niccolai2020/Niccolai2020_output_from_MSI/")
#NOTE: change this when actual files come in. 

filelist = glob.glob('*_V01_Niccolai_input_output.csv') #why haven't I been doing this forever? this creates a list of files with the specified string

#make a list of sample names, a list of dataframes, and a dictionary of both of these together. 
sample_names_list = []
for file in filelist:
    sample_name = re.sub("_V01_Niccolai_input_output.csv", "", file) #delete extra crap on string.  
    sample_names_list.append(sample_name)
del file, sample_name

df_list = []
for file in filelist:
    df= pd.read_csv(file)
    df_list.append(df)
del file, df

df_dict = dict(zip(sample_names_list, df_list))

#####step 2: get a list of all unique metabolite values in each file
all_metab_values = []
for df in df_list:
    metab_list = df["Unnamed: 0"].tolist()
    all_metab_values.extend(metab_list)
del df, metab_list

unique_metab_values = set(all_metab_values) 

#cleanup point
del all_metab_values, df_list, filelist, sample_names_list

#clean up metab values (remove terrible prefixes and suffixes)
unique_metabs_shortname_list = []
for item in unique_metab_values:
    item= item.replace("EX_", "")
    item= item.replace("_m", "")
    unique_metabs_shortname_list.append(item)
del item

unique_metabs_dict = dict(zip(unique_metab_values, unique_metabs_shortname_list))

#####step 3: create a dataframe to append other dataframes to.
combo_metab_df = pd.DataFrame.from_dict(unique_metabs_dict, orient='index')
combo_metab_df.reset_index(inplace=True)
combo_metab_df.rename(columns={"index": "model_metabolite_names", 0: "abbreviated_metabolite_names"}, inplace=True)
#note here that we need both model_metabolite_names (contains EX_ and _m) and abbreviated_metabolite_names (ends trimmed)
#the former is for merging in data; the latter is for assigning longform names

#cleanup
del unique_metabs_dict

#####step 4: append all input_output files to combo_df
for key, value in df_dict.items():
    value.rename(columns={"Unnamed: 0": "metabolite_names", "0": "metabolite_values"}, inplace=True)
    indiv_metabs_dict= value.set_index('metabolite_names').to_dict()['metabolite_values']
    combo_metab_df[str(key)] = combo_metab_df['model_metabolite_names'].map(indiv_metabs_dict)
del indiv_metabs_dict, key, value    

#####step 5: rename metabolites into something useful.
#from 10.12.20_community_medium_data_processing.py
#step 5a: import metabolites key. NEW KEY
#hopefully the Burns one works fine for Hale data also.
metab_key = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Common_growth_medium_runs/May_June2021_MICOM_flux_and_pathway_analyses/Metabolites_key_full_May2021_updates_V2.csv", sep=',')    #read in the metabolite names key
#there are a lot of extra columns here: ony really want the first 2
metab_key_slimmed = metab_key.loc[:, 'abbreviation':'fullName']    #use : to select full column
metab_key_slimmed = metab_key_slimmed.rename(columns={"abbreviation": "Abbreviation", "fullName": "FullName"})
del metab_key
  
#step 5b: create a dictionary with unique_metab_values and their associated long names
#look up all values of unique_metab_values in metab_key_slimmed and create a dictionary
unique_metabs_long_plus_short = {}

metab_key_slimmed_dict= metab_key_slimmed.set_index('Abbreviation').to_dict()['FullName'] #convert metab_key_slimmed to dict for lookup
for key, value in metab_key_slimmed_dict.items():
    if str(key) in unique_metabs_shortname_list:
        #print(key)
        unique_metabs_long_plus_short[key] = value
del key, value

#step 5c: append the shortform metab names to the combo_metab_df as a column using our old friend map.
combo_metab_df["longform_metabolite_names"] = combo_metab_df['abbreviated_metabolite_names'].map(unique_metabs_long_plus_short)

#cleanup
del df_dict, unique_metab_values, unique_metabs_long_plus_short, unique_metabs_shortname_list

#double check here that all metabolites have their names assigned. (i.e. all longform_metabolite_name row values are present)
#Looks like we're all good here! Include extra steps but comment them out.

#############################Missing metabolite discovery INTERM START################################
#goal: update Metabolites_key_full.tsv with metabolites not present in this list but identified in the models
#nan values in longform_metabolite_names are the ones to look up. 
#nan values in sample columns are metabolites not imported/exported in that community
#os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Hale2018/")

#combo_metab_df_for_lookup = combo_metab_df[combo_metab_df["longform_metabolite_names"].isnull()]
#combo_metab_df_for_lookup = combo_metab_df_for_lookup.loc[:, ["model_metabolite_names", "abbreviated_metabolite_names", "longform_metabolite_names"]]
#combo_metab_df_for_lookup.to_csv("09.02.21_unknown_Burns_metabolites_for_lookup.csv")
#Look these metabolites up, hand add them, and re-import the resulting dataframe
#missing_metabs= pd.read_csv("09.02.21_unknown_Burns_metabolites_for_lookup_lookupdone.csv")
#missing_metabs_dict = dict(zip(missing_metabs.abbreviated_metabolite_names, missing_metabs.longform_metabolite_names))
#NEW: add missing_metabs_dict to metab_key_slimmed- need this later
#metab_key_missing= pd.DataFrame.from_dict(missing_metabs_dict, orient='index')
#metab_key_missing.reset_index(inplace=True)
#metab_key_missing= metab_key_missing.rename(columns={"index": "Abbreviation", 0: "FullName"})
#metab_key_slimmed_plus_missing = pd.concat([metab_key_slimmed, metab_key_missing], axis=0)
#now back to old stuff. 
#merge this df into combo_metab_df
#combo_metab_df2= copy.deepcopy(combo_metab_df)
#combo_metab_df2["longform_metabolite_names2"] = combo_metab_df2['abbreviated_metabolite_names'].map(missing_metabs_dict)
#now we have two columns for longform metabolite names. Need to merge them without losing data
#combo_metab_df2["longform_metabolite_names"] = np.where(combo_metab_df2["longform_metabolite_names"].isna(),\
                                                        #combo_metab_df2["longform_metabolite_names2"], \
                                                            #combo_metab_df2["longform_metabolite_names"]).astype("str")
#basically, look through col1 and fill in the values from col2 wherever col1 is nan.
#and, of course, delete col2    
#combo_metab_df2 = combo_metab_df2.drop("longform_metabolite_names2", axis=1)
#############################Missing metabolite discovery INTERM OVER################################

#step 6: clean up the dataframe and prepare for export
combo_metab_df3= copy.deepcopy(combo_metab_df)
combo_metab_df3= combo_metab_df3.set_index('longform_metabolite_names')
combo_metab_df3.drop(columns=["model_metabolite_names", "abbreviated_metabolite_names"], inplace=True)
combo_metab_df3.reset_index(inplace=True) #reset index (temporarily)
combo_metab_df3= combo_metab_df3.set_index('longform_metabolite_names')
#combo_metab_df3.to_csv("09.02.21_Hale2018_metabolites_imported_exported.csv")

#step 7: traspose and add metadata.
#first, transpose
combo_metab_df3_flipped = combo_metab_df3.T #set index and tranpose.
combo_metab_df3_flipped.reset_index(inplace=True)
combo_metab_df3_flipped_longform = combo_metab_df3_flipped.melt(id_vars=["index"], var_name="metabolite_name", value_name="metabolite_amount")
#here, "index" is columns to keep (id_vars), "fluxpath_name" is category to merge under (var_name),
#and "fluxpath_amount" is name of column with category values (value_name)    
combo_metab_df3_flipped_longform = combo_metab_df3_flipped_longform.rename(columns={"index": "SampleID"}) #name sample ID column
#then, add metabolite metadata
combo_metab_df3_flipped_longform = pd.merge(left=combo_metab_df3_flipped_longform , right=metab_key_slimmed, how='left', 
                                          left_on='metabolite_name', right_on='FullName')
#if doing the interm for metabolite discovery, right merge is on metab_key_slimmed_plus_missing
combo_metab_df3_flipped_longform = combo_metab_df3_flipped_longform.drop(columns=["Abbreviation", "FullName"])

#then, add sample metadata
metadata = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Niccolai2020/Niccolai2020_all_metadata.csv")
#wtf, why are therea bunch of extra rows here?
#metadata.dropna(subset = ["Patient_Blind_ID"], inplace=True)
combo_metab_df3_plus_metadata = pd.merge(left=combo_metab_df3_flipped_longform, right=metadata, 
                                                        how='left', left_on='SampleID', right_on='Sample_ID') #merge metadata
combo_metab_df3_plus_metadata.drop(['SampleID', 'Age', 'Diagnosis', 'TNM', 'Statium', 'Site'], axis=1, inplace=True)#drop superfluous columns
#export!
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Niccolai2020/") #analysis directory
combo_metab_df3_plus_metadata.to_csv("11.09.21_Niccolai2020_inputoutput_with_metadata_longform.csv")











































