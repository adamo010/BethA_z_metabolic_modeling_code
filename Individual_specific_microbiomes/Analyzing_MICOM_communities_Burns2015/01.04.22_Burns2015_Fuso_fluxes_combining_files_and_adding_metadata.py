#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 16:46:07 2022

@author: adamo010
"""
##newerer starting point: Date_TBA_combining_flux_values_from_MSI_runs.py; adapting this from Burns to Hale data
#(and from 09.02.21_Hale2018_fluxes_combining_files_and_adding_metadata.py)
#newer starting point: 06.16.21_combinding_TDpointonefive_flux_values_from_MSI_runs.py
#starting point: 05.26.21_combining_TDpointone_inputoutput_values_from_MSI_runs.py

#goal: take all fluxes from MSI run, and organize them into something that can be correlated.
#update from 09.09.21: ONLY using fuso fluxes. Stupid.f

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

#######step 2: find Fuso OTU_IDs
#first, need to find the fuso OTU_IDs
filtered_tax = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/09.09.21_Burns2015_indiv_spp_GRs_collapsed_to_genus.csv")
fuso_only_tax =  filtered_tax[filtered_tax['genus']=="Fusobacterium"]
fuso_OTU_IDs = fuso_only_tax['OTU_ID'].unique() #get unique OTU_IDs of fuso.
for elem in fuso_OTU_IDs:
    fuso_ID = str(elem)
del(fuso_OTU_IDs, elem)

########step 3: transpose and filter by fuso-only taxa
df_dict = dict(zip(sample_names_list, df_list)) #add back later, for transposed DFs

df_list2 = []
for key, value in df_dict.items():
    df=value.set_index('compartment')
    df= df.drop("medium", axis=0)
    df_transposed = df.transpose()
    df_transposed['SampleID']=str(key)
    df_transposed.drop(df_transposed.columns.difference([fuso_ID, "SampleID"]),1,inplace=True)
    df_list2.append(df_transposed)
del df, key, value, df_transposed    
    
########step 4: merge all dataframes and filter out fuso-negative sample IDs
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

#########step 5: filter out fluxpaths which are 0 or NAN in ALL Fuso-positive samples
#aw cripes, how do I do this? probably going to have to count all the nans or something.
#all right, hang on. Let's get a list of all the flux names
combo_df_filtered.reset_index(inplace=True)
metabs_list = combo_df_filtered['index'].unique().tolist() 
combo_df_filtered.rename(columns={"index":"flux_name"}, inplace=True)

#count the number of nan values for each fluxpath name
metabs_nan_count = []
for item in metabs_list:
    df1 = combo_df_filtered[combo_df_filtered['flux_name'] == item] #create a dataframe subset containing a single flux_name
    count = df1['flux_value'].isnull().sum()
    metabs_nan_count.append(count)
del item, df1, count    

#create a dictionary of fluxpaths and their nan counts
metabs_dict = dict(zip(metabs_list, metabs_nan_count))    

#create two dictionaries, one which contains all fluxpaths with nans in all 44 samples (all fuso positive samples), and the rest
metabs_nan_dict = {}
metabs_nonnan_dict = {}
for key, value in metabs_dict.items():
    if value == 44:
        metabs_nan_dict[key]=value
    elif value <= 43:
        metabs_nonnan_dict[key]=value
del key, value

#now, create a list of all fluxpaths all fluxpaths with nans in all 44 samples and use it to filter combo_df_filtered
nonnan_metabs_list = list(metabs_nonnan_dict.keys())
combo_df_filtered2 = combo_df_filtered[combo_df_filtered['flux_name'].isin(nonnan_metabs_list)] 
     
######step 6: rename fluxes into something useful.
#from 10.12.20_community_medium_data_processing.py
#step 5a: import metabolites key. NEW KEY
fluxpath_key = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Common_growth_medium_runs/May_June2021_MICOM_flux_and_pathway_analyses/Microbiome_flux_key_full_June2021.csv", sep=',')#there are a lot of extra columns here: ony really want the first 2
fluxpath_key_slimmed = fluxpath_key.loc[:, 'abbreviation':'subsystem']    #use : to select full column
del fluxpath_key

#combo_df_filtered.reset_index(inplace=True) #reset index
combo_df_filtered2.rename(columns={"flux_name": "fluxpath_name"}, inplace=True) #name column with fluxpath ids as full_fluxpath_names
combo_df_plus_metab_metadata = pd.merge(left=combo_df_filtered2, right=fluxpath_key_slimmed, how='left', 
                                          left_on='fluxpath_name', right_on='abbreviation')
#hmmmmmm, that worked, but there are a lot of missing data here. I wonder if it's because we did exactly zero filtering.
#do we need to filter? I'm not sure. We do catch the zero values in the stats file in R, but...
#I think I'll try to add filtering back in.

#then, add sample metadata
metadata = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Burns_2015_sample_metadata.csv")#wtf, why are therea bunch of extra rows here?
metadata.dropna(subset = ["Patient_Blind_ID"], inplace=True)
combo_fluxpath_df_plus_metab_sample_metadata = pd.merge(left=combo_df_plus_metab_metadata, right=metadata, 
                                                        how='left', left_on='SampleID', right_on='SampleID') #merge metadata
combo_fluxpath_df_plus_metab_sample_metadata.drop(['Site', 'MSI_status', 'Stage'], axis=1, inplace=True)#drop superfluous columns
#export!
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/") #analysis directory
combo_fluxpath_df_plus_metab_sample_metadata.to_csv("01.04.22_Burns2015_Fuso_fluxes_with_metadata_longform.csv")


#######################JUNK/testing

#
testo = df_list[1]
testo=testo.set_index('compartment').drop("medium", axis=0)
#
testo= testo[testo['compartment']== str(fuso_ID)]
testo=testo.set_index('compartment')
testo=testo.drop("medium", axis=0)
testo_transposed = testo.transpose()
testo_transposed["sample_ID"]=str("testid")
testo_transposed.drop(testo_transposed.columns.difference([fuso_ID, "sample_ID"]),1,inplace=True)







#then, filter each df such that only fuso OTU_ids are included.
df_list2 = []
for df in df_list:
    df2 = df[df['compartment'] == str(fuso_ID)]
    df_list2.append(df2)
del df, df2
#note! now, there will be several empty dfs, as these samples lack fuso. How... will we account for this? worry about it later. 

#do I need this?- how to get Fuso-containing PBIDs. Shouldn't need this, if we remove all non-fuso OTU-ids, but might be a good sanity check.
fuso_pos_samples = fuso_only_tax.loc[:, ['Patient_Blind_ID', 'SampleID']] #only keep two relevant columns
paired_fuso_PBIDs = fuso_pos_samples[fuso_pos_samples.duplicated('Patient_Blind_ID', keep=False) == True] #only keep PBIDs with paired samples
fuso_PBIDs_prelist = paired_fuso_PBIDs['Patient_Blind_ID'].astype(str).unique().tolist() #convert PBIDs to list
fuso_PBIDs_list = []
for item in fuso_PBIDs_prelist:
    item2 = item[:-2]
    fuso_PBIDs_list.append(item2)
del(item, item2, filtered_tax, fuso_only_tax, fuso_pos_samples, fuso_PBIDs_prelist, paired_fuso_PBIDs)

########step 3: stack all the dataframes together.
combo_df = []
for df in df_list:
   







