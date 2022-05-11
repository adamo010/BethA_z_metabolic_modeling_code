#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 13:45:20 2022

@author: adamo010
"""

##newerer starting point: Date_TBA_combining_flux_values_from_MSI_runs.py; adapting this from Burns to Hale data
#newer starting point: 06.16.21_combinding_TDpointonefive_flux_values_from_MSI_runs.py
#starting point: 05.26.21_combining_TDpointone_inputoutput_values_from_MSI_runs.py
#update from 09.02.21 version: only moving forward with EX_ metabolites.
#borrowed heavily from 10.13.21_Hale2018_fluxes_combining_files_and_adding_metadata.py


#goal: take all fluxes from MSI run, and organize them into something that can be correlated.
#update from 11.09.21: ONLY using fuso fluxes. Stupid.f

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
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Niccolai2020/Niccolai2020_output_from_MSI/")
#NOTE: change this when actual files come in. 

filelist = glob.glob('*_V01_Niccolai_fluxes.csv') #why haven't I been doing this forever? this creates a list of files with the specified string

sample_names_list = []
for file in filelist:
    sample_name = re.sub("_V01_Niccolai_fluxes.csv", "", file) #delete extra crap on string.  
    sample_names_list.append(sample_name)
del file, sample_name

df_list = []
for file in filelist:
    df= pd.read_csv(file)
    df_list.append(df)
del file, df

df_dict = dict(zip(sample_names_list, df_list)) #add back later, for transposed DFs

#testing areas
#test_sample_name= sample_name
#test_df = df

######step 2: clean up files and prepare for merging
df_list_cleaned = []
for key, value in df_dict.items():
    df_cleaned = value.rename(columns={"compartment": "OTU_ID"})
    df_cleaned = df_cleaned.set_index('OTU_ID')
    df_cleaned= df_cleaned.drop("medium", axis=0)
    df_transposed = df_cleaned.T
    df_transposed.reset_index(inplace=True)
    df_transposed = df_transposed.rename(columns={"index": "fluxpath_id"})
    df_long =pd.melt(df_transposed,id_vars=['fluxpath_id'],var_name='otu_id', value_name='fluxpath_values')
    df_long['SampleID']=str(key)
    df_list_cleaned.append(df_long)
del key, value, df_cleaned, df_transposed, df_long    

#######step 2: find Fuso OTU_IDs
#first, need to find the fuso OTU_IDs
filtered_tax = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Niccolai2020/Niccolai_2020_OTUs_with_matched_taxonomy_V3.csv")
fuso_only_tax =  filtered_tax[filtered_tax['Genus']=="Fusobacterium"]
fuso_OTU_IDs = fuso_only_tax['OTU_ID'].unique() #get unique OTU_IDs of fuso.
fuso_IDs_list=[]
for elem in fuso_OTU_IDs:
    #fuso_ID = str(elem)
    fuso_IDs_list.append(elem)
del(fuso_OTU_IDs, elem)

sorting_df = filtered_tax[["Genus", "OTU_ID"]].copy()

#####step 3 filter by fuso-only otu ids

#full code
df_list_filtered = []
for item in df_list_cleaned:
    df_filtered= pd.merge(left=item, right= sorting_df, how='left', left_on='otu_id', right_on="OTU_ID")
    df_filtered = df_filtered[df_filtered.Genus == "Fusobacterium"]
    df_list_filtered.append(df_filtered)
del(item, df_filtered) 

#then stak. all column names should be the same
########step 4: merge all dataframes
combo_df = pd.concat(df_list_filtered, axis=0)
#nice
combo_df = combo_df.drop(columns={"otu_id"}) #drop extra column

########step 5: Add PBIDs and remove unpaired samples
metadata = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Niccolai2020/Niccolai2020_all_metadata.csv")
#metadata.dropna(axis=0, thresh=13, inplace=True) #remove blank rows; thresh=13 b/c there are 13 columns, and we only want to remove rows with nan in all columns
rel_metadata = metadata[["Sample_ID","Patient_Blind_ID","Description"]]
#rel_metadata.rename(columns={"sampleIDs_from_read_counts": "SampleID", "host_subject_id":"Patient_Blind_ID",\
                             #"normal_adjacent_or_tumor_tissue_specimen": "Description"}, inplace=True)

#which sampleids are fuso positive?
filtsampleIDlist = combo_df['SampleID'].astype(str).unique().tolist() #55 samples have Fuso. 

pos_list = []
for row in rel_metadata["Sample_ID"]:
    if row in filtsampleIDlist:
        Fuso = "yes"
    else:
        Fuso="no"
    pos_list.append(Fuso)
del row, Fuso
rel_metadata['Fuso_pos']= pos_list #append pos-list to metadata frame
del pos_list

paired_rel_metadata = rel_metadata[rel_metadata.duplicated(subset=['Patient_Blind_ID','Fuso_pos'], keep=False) == True] #remove discordant samples
paired_rel_metadata = paired_rel_metadata[paired_rel_metadata.Fuso_pos == "yes"] #remove PBIDs where neither sample has fuso
paired_sampleids = paired_rel_metadata["Sample_ID"].unique().tolist() #convert SampleIDs to list. There are 72. That is correct? Maybe not. But close enough

#now, we filter combo_df by paired_sampleids
combo_df_fuso = combo_df[combo_df["SampleID"].isin(paired_sampleids)]
combo_fuso_df_with_pbid = pd.merge(left=combo_df_fuso, right= paired_rel_metadata, how="left", left_on="SampleID", right_on="Sample_ID")
combo_fuso_df_with_pbid.drop(columns={"Fuso_pos"}, inplace=True)
#good

######step 6: rename fluxes into something useful.
#from 10.12.20_community_medium_data_processing.py
#step 5a: import metabolites key. NEW KEY
fluxpath_key = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Common_growth_medium_runs/May_June2021_MICOM_flux_and_pathway_analyses/Microbiome_flux_key_full_June2021.csv", sep=',')#there are a lot of extra columns here: ony really want the first 2
fluxpath_key_slimmed = fluxpath_key.loc[:, 'abbreviation':'subsystem']    #use : to select full column
del fluxpath_key

combo_df_plus_metab_metadata = pd.merge(left=combo_fuso_df_with_pbid, right=fluxpath_key_slimmed, how='left', 
                                          left_on='fluxpath_id', right_on='abbreviation')

#export!
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Niccolai2020/") #analysis directory
combo_df_plus_metab_metadata.to_csv("01.10.22_Niccolai2020_Fuso_fluxes_with_metadata_longform.csv")
























