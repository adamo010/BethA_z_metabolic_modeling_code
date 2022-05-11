#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 13:03:32 2022

@author: adamo010
"""

##newerer starting point: Date_TBA_combining_flux_values_from_MSI_runs.py; adapting this from Burns to Hale data
#newer starting point: 06.16.21_combinding_TDpointonefive_flux_values_from_MSI_runs.py
#starting point: 05.26.21_combining_TDpointone_inputoutput_values_from_MSI_runs.py
#update from 09.02.21 version: only moving forward with EX_ metabolites.


#goal: take all fluxes from MSI run, and organize them into something that can be correlated.
#update from 09.02.21: ONLY using fuso fluxes. Stupid.

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

######step 1: import filesle2018_MICOM_output_from_MSI/Hale_V12/")
#NOTE: change this when actual files come in. 
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Running_MICOM_communities_Hale2018/Hale2018_MICOM_output_from_MSI/Hale_V12/")
filelist = glob.glob('*_V12_Hale_fluxes.csv') 

sample_names_list = []
for file in filelist:
    sample_name = re.sub("_V12_Hale_fluxes.csv", "", file) #delete extra crap on string.  
    sample_names_list.append(sample_name)
del file, sample_name

df_list = []
for file in filelist:
    df= pd.read_csv(file)
    df_list.append(df)
del file, df

df_dict = dict(zip(sample_names_list, df_list)) #add back later, for transposed DFs

#testing areas
test_sample_name= sample_name
test_df = df

######step 2: clean up files and prepare for merging

#testing areas
test_df_cleaned = test_df.rename(columns={"compartment": "OTU_ID"})
test_df_cleaned=test_df_cleaned.set_index('OTU_ID')
test_df_cleaned= test_df_cleaned.drop("medium", axis=0)
test_df_transposed = test_df_cleaned.T
test_df_transposed.reset_index(inplace=True)
test_df_transposed = test_df_transposed.rename(columns={"index": "fluxpath_id"})
test_df_long =pd.melt(test_df_transposed,id_vars=['fluxpath_id'],var_name='otu_id', value_name='fluxpath_values')
test_df_long['SampleID']=str(test_sample_name)

#full code
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
filtered_tax = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Hale2018/Hale_2018_OTUs_with_matched_taxonomy_and_proxy_OTUids.csv")
fuso_only_tax =  filtered_tax[filtered_tax['Genus']=="Fusobacterium"]
fuso_OTU_IDs = fuso_only_tax['OTU_ID'].unique() #get unique OTU_IDs of fuso.
fuso_IDs_list=[]
for elem in fuso_OTU_IDs:
    #fuso_ID = str(elem)
    fuso_IDs_list.append(elem)
del(fuso_OTU_IDs, elem)

sorting_df = filtered_tax[["Genus", "OTU_ID"]].copy()

#####step 3 filter by fuso-only otu ids

#testing areas
test_df_filtered= pd.merge(left=test_df_long, right= sorting_df, how='left', left_on='otu_id', right_on="OTU_ID")
test_df_filtered = test_df_filtered[test_df_filtered.Genus == "Fusobacterium"]

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
metadata = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Hale2018/Hale_2018_metadata_for_MICOM_samples.csv")
metadata.dropna(axis=0, thresh=13, inplace=True) #remove blank rows; thresh=13 b/c there are 13 columns, and we only want to remove rows with nan in all columns
rel_metadata = metadata[["sampleIDs_from_read_counts","host_subject_id","normal_adjacent_or_tumor_tissue_specimen"]]
rel_metadata.rename(columns={"sampleIDs_from_read_counts": "SampleID", "host_subject_id":"Patient_Blind_ID",\
                             "normal_adjacent_or_tumor_tissue_specimen": "Description"}, inplace=True)

#which sampleids are fuso positive?
filtsampleIDlist = combo_df['SampleID'].astype(str).unique().tolist() #55 samples have Fuso. 

pos_list = []
for row in rel_metadata["SampleID"]:
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
paired_sampleids = paired_rel_metadata["SampleID"].unique().tolist() #convert SampleIDs to list. There are 46. That is correct. 

#now, we filter combo_df by paired_sampleids
combo_df_fuso = combo_df[combo_df["SampleID"].isin(paired_sampleids)]
combo_fuso_df_with_pbid = pd.merge(left=combo_df_fuso, right= paired_rel_metadata, how="left", left_on="SampleID", right_on="SampleID")
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
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Hale2018/") #analysis directory
combo_df_plus_metab_metadata.to_csv("01.07.22_Hale2018_Fuso_fluxes_with_metadata_longform.csv")



    






##########################JUNK###############

    
###########add Genus onto each df and filter by fuso only taxa

df_list2 = []
for key, value in df_dict.items():
    app_df = pd.merge(left=value, right= sorting_df, how='left', left_on='compartment', right_on="OTU_ID")
    app_df['SampleID']=str(key)
    filtered_df = app_df[app_df.Genus == "Fusobacterium"]
    df_list2.append(filtered_df)
del(app_df, filtered_df, key, value)
    
########step 3: transpose
df_list3 = []

for item in df_list:
    df=item.set_index('compartment')
    #df= df.drop("medium", axis=0)
    df_transposed = df.transpose()
    #df_transposed['SampleID']=str(key)
    #df_transposed.drop(df_transposed.columns.difference([fuso_ID, "SampleID"]),1,inplace=True)
    df_transposed.reset_index(inplace=True)
    df_list3.append(df_transposed)
#del df, item, df_transposed    

test_sample_df = df_transposed  

#####NEW step 4/old step 3: get a list of all unique flux pathway values in each file
all_fluxpath_values = []
for df in df_list:
    ph = df
    ph=ph.drop(columns="compartment")
    #df.set_index("compartment",inplace=True)
    df_transposed = ph.T
    df_transposed.reset_index(inplace=True)
    fluxpath_list = df_transposed["index"].tolist()
    all_fluxpath_values.extend(fluxpath_list)
del df, fluxpath_list

unique_fluxpath_values = set(all_fluxpath_values) 
#INTERESTING:  192 fluxpaths for Burns, 324 for Hale

#CLEANUP
del all_fluxpath_values, df_list, filelist, sample_names_list

######NEW step 5/ old step 4: create a dataframe to append other dataframes to.
#don't actually need this step here.
combo_fluxpath_df = pd.DataFrame(unique_fluxpath_values)
combo_fluxpath_df.rename(columns={0: "abbreviated_fluxpath_names"}, inplace=True)
#note here that we need both model_metabolite_names (contains EX_ and _m) and abbreviated_metabolite_names (ends trimmed)
#the former is for merging in data; the latter is for assigning longform names

for item in df_list3:
    combo_fluxpath_df= pd.merge(left=combo_fluxpath_df, right=item, how="left", left_on="abbreviated_fluxpath_names", right_on="index")

test_merged_df = pd.merge(left=combo_fluxpath_df, right= test_sample_df, how="left", left_on= "abbreviated_fluxpath_names", right_on="index")

#this didn't work. Start here. 

########step 4: merge all dataframes
combo_df = pd.concat(df_list2, axis=0)
combo_df.rename(columns={"1096339":"flux_value"}, inplace=True)
#how many sample IDs? 
sampleIDlist = combo_df['SampleID'].astype(str).unique().tolist() #88 samples, so all samples included here (which was the goal)
#now, we get a list of Fuso-containing sampleIDs
fuso_pos_samples = fuso_only_tax.loc[:, ['Patient_Blind_ID', 'SampleID']] #only keep two relevant columns
paired_fuso_PBIDs = fuso_pos_samples[fuso_pos_samples.duplicated('Patient_Blind_ID', keep=False) == True] #only keep PBIDs with paired samples
fuso_sampleids_list = paired_fuso_PBIDs['SampleID'].astype(str).unique().tolist() #convert SampleIDs to list
#if needed, to get paired fuso PBIDs
#fuso_PBIDs_prelist = paired_fuso_PBIDs['Patient_Blind_ID'].astype(str).unique().tolist() #convert PBIDs to list
#fuso_PBIDs_list = []
#for item in fuso_PBIDs_prelist:
    #item2 = item[:-2]
    #fuso_PBIDs_list.append(item2)
#del(item, item2, filtered_tax, fuso_only_tax, fuso_pos_samples, fuso_PBIDs_prelist, paired_fuso_PBIDs)
#now, filter combo_df by fuso_sampleids_list
boolean_series = combo_df.SampleID.isin(fuso_sampleids_list)
combo_df_filtered = combo_df[boolean_series]
del boolean_series
#check again- how many unique Sampleids do we have? should match the number of sampleids in fuso_sampleids_list (here, 44)
filtsampleIDlist = combo_df_filtered['SampleID'].astype(str).unique().tolist() #88 samples, so all samples included here (which was the goal)
#nice
del sampleIDlist, filtsampleIDlist















    
