#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 12:34:03 2021

@author: adamo010
"""

##newerer starting point: Date_TBA_combining_flux_values_from_MSI_runs.py; adapting this from Burns to Hale data
#newer starting point: 06.16.21_combinding_TDpointonefive_flux_values_from_MSI_runs.py
#starting point: 05.26.21_combining_TDpointone_inputoutput_values_from_MSI_runs.py
#update from 09.02.21 version: only moving forward with EX_ metabolites.


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
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Running_MICOM_communities_Hale2018/Hale2018_MICOM_output_from_MSI/Hale_V12/")
#NOTE: change this when actual files come in. 

filelist = glob.glob('*_V12_Hale_fluxes.csv') #why haven't I been doing this forever? this creates a list of files with the specified string

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

####NEW NEW step 2: transpose and only keep rows with "EX_" in their names (exchange Rxs)
df_list2a = []

for df in df_list:
    df=df.set_index('compartment') #set row names to flux pathway names  
    df=df.drop("medium", axis=0) #drop the medium column (not exactly sure what it means in terms of flux anyway)
    df_transposed = df.transpose() #flipperoo so flux pathway names are row names, compartments (taxa) are column names
    df_transposed.reset_index(inplace=True) #reste index
    df_transposed= df_transposed.rename(columns={"index": "fluxpath"}) #rename index column as fluxpath
    df_transposed2 = df_transposed[df_transposed.fluxpath.str.contains("EX_", na=False)] #only keep fluxpaths with EX_ in name; need na=False to run with NA values in df
    df_transposed2.set_index("fluxpath", inplace=True)
    otu_names = []
    for col in df_transposed2.columns:       #get a list of otu_names
        otu_names.append(col)
    del col
    num_above_cutoff = [(abs(df_transposed2[otu_names]) >= 0.000001).sum(1)] #count the number of OTUs with flux >=10-6 for each fluxpath
    #convert to a list
    for elem in num_above_cutoff:
        otu_table_num_above_cutoff = elem.values.tolist()       #this converts a series to a list
    del elem
    df_transposed2['abs_Num_above_cutoff'] = otu_table_num_above_cutoff #add this list as a column
    num_otus = len(otu_names) #get number of OTUs in this sample
    df_transposed2_above_cutoff = df_transposed2.loc[df_transposed2['abs_Num_above_cutoff'] >= (num_otus*0.1)]
    #only keep fluxpaths that appear in at least 10% of OTUs at a flux value of >=10-6
    df_transposed2_above_cutoff.reset_index(inplace=True) #reset index
    df_transposed2_above_cutoff= df_transposed2_above_cutoff.drop(["abs_Num_above_cutoff"], axis=1) #drop cutoff column
    df_long = df_transposed2_above_cutoff.melt(id_vars= ["fluxpath"], var_name = "OTU_ID", value_name= "flux_value")
    #keeping fluxpath column constant, merge all other columns into 2 columns: OTU_id (old col names) and flux values (old data values)
    df_list2a.append(df_long)
    del df_long, df_transposed, df_transposed2, df_transposed2_above_cutoff, num_above_cutoff, num_otus, otu_names, otu_table_num_above_cutoff
#all right, that took for fucking ever. Moving on.
    
df_dict = dict(zip(sample_names_list, df_list2a)) #add back later, for transposed DFs

####NEW step 3: add sample_names to each DF for later merging.
for key, value in df_dict.items():
    sample_id= str(key)
    value['sample_name'] = pd.Series([sample_id for x in range(len(value.index))]) #add the sample ID as a constant in a new col in each df
del key, value, sample_id
#neat.

#####NEW step 4/old step 3: get a list of all unique flux pathway values in each file
testo= df_list2a[0]

all_fluxpath_values = []
for df in df_list2a:
    #df = df.reset_index()
    fluxpath_list = df["fluxpath"].tolist()
    all_fluxpath_values.extend(fluxpath_list)
del df, fluxpath_list

unique_fluxpath_values = set(all_fluxpath_values) 
#INTERESTING:  192 fluxpaths for Burns, 324 for Hale

#CLEANUP
del all_fluxpath_values, df_list, filelist, sample_names_list

######NEW step 5/ old step 4: create a dataframe to append other dataframes to.
#don't actually need this step here.
#combo_fluxpath_df = pd.DataFrame(unique_fluxpath_values)
#combo_fluxpath_df.reset_index(inplace=True)
#combo_fluxpath_df.rename(columns={0: "abbreviated_fluxpath_names"}, inplace=True)
#note here that we need both model_metabolite_names (contains EX_ and _m) and abbreviated_metabolite_names (ends trimmed)
#the former is for merging in data; the latter is for assigning longform names

######NEW step 6/old step 5: append all input_output files to combo_df
# a change from the old version here: I think we can just stack the dataframes. LOOOOOOOONG version.
dfs = df_dict.values()
df_list3= list(dfs)
combo_fluxpath_df = pd.concat(df_list3)

#merged_left = pd.merge(left=survey_sub, right=species_sub, how='left', left_on='species_id', right_on='species_id')
#result = pd.concat(frames)


######step 6: rename fluxes into something useful.
#from 10.12.20_community_medium_data_processing.py
#step 5a: import metabolites key. NEW KEY
fluxpath_key = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Common_growth_medium_runs/May_June2021_MICOM_flux_and_pathway_analyses/Microbiome_flux_key_full_June2021.csv", sep=',')#there are a lot of extra columns here: ony really want the first 2
fluxpath_key_slimmed = fluxpath_key.loc[:, 'abbreviation':'subsystem']    #use : to select full column
del fluxpath_key

######step 7: transpose and add sample metadata.
#first, transpose and convert to longform. ALREADY DONE.
#combo_fluxpath_df_flipped = combo_fluxpath_df.set_index("abbreviated_fluxpath_names").T #set index and tranpose.
#combo_fluxpath_df_flipped.reset_index(inplace=True)
#combo_fluxpath_df_flipped_longform = combo_fluxpath_df_flipped.melt(id_vars=["index"], var_name="fluxpath_name", value_name="fluxpath_amount")
#here, "index" is columns to keep (id_vars), "fluxpath_name" is category to merge under (var_name),
#and "fluxpath_amount" is name of column with category values (value_name)    
combo_fluxpath_df = combo_fluxpath_df.rename(columns={"sample_name": "Sample_ID"}) #name sample ID column
#then, add metabolite metadata
combo_fluxpath_plus_metab_metadata = pd.merge(left=combo_fluxpath_df, right=fluxpath_key_slimmed, how='left', 
                                          left_on='fluxpath', right_on='abbreviation')
combo_fluxpath_plus_metab_metadata = combo_fluxpath_plus_metab_metadata.drop(columns=["abbreviation"]) #redundant
combo_fluxpath_plus_metab_metadata = combo_fluxpath_plus_metab_metadata.rename(columns={"description": "fluxpath_description",
                                                                                        "formula": "fluxpath_formula",
                                                                                        "subsystem": "fluxpath_subsystem",
                                                                                        "fluxpath": "fluxpath_name"}) #name sample ID column

#then, add sample metadata
metadata = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Hale2018/Hale_2018_metadata_for_MICOM_samples.csv")
metadata.dropna(subset = ["sampleIDs_from_read_counts"], inplace=True)
combo_fluxpath_df_plus_metab_sample_metadata = pd.merge(left=combo_fluxpath_plus_metab_metadata, right=metadata, 
                                                        how='left', left_on='Sample_ID', right_on='sampleIDs_from_read_counts') #merge metadata
combo_fluxpath_df_plus_metab_sample_metadata.drop(["cancer_stage", "env_material", "ethnicity", "Host_Age", "host_body_mass_index",\
           "host_sex", "host_tissue_sampled", "mismatch_repair_status", "collection_site", "biopsy_number"], axis=1, inplace=True)#drop superfluous columns
#export!
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Hale2018/")
combo_fluxpath_df_plus_metab_sample_metadata.to_csv("10.13.21_Hale2018_fluxes_with_metadata_longform.csv")


