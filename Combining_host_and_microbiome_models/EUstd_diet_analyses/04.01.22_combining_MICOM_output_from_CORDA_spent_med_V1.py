#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 12:37:54 2022

@author: adamo010
"""

import scipy
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

#the goal of this is to combine a bunch of output files from CORDA spent medium --> MICOM models
#we want to combine four different types of files:
#Part 1: CORDA spent medium
#Part 2: MICOM community growth rates:
#Part 3: MICOM indiv spp growth rates: see 11.05.21_Niccolai2020_indiv_spp_GRs_combining_files_and_adding_metadata.py
#Part 4: MICOM fluxes

############Part 1: CORDA spent medium
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_host_spent_med_localgen/")
spentmed_file_list = glob.glob('*hostmodel_EUstd_spent_medium.csv')

spentmed_file_names = []
for elem in spentmed_file_list:
    sample_name = re.sub("_hostmodel_EUstd_spent_medium.csv", "", elem)  
    spentmed_file_names.append(sample_name)        
del elem, sample_name    

df_list = []
for file in spentmed_file_list:
    df= pd.read_csv(file)
    df=df.drop(columns={"Unnamed: 0"})
    df=df.rename(columns={"premed_name": "metabolite_name", "net_postmed_amount": "metabolite_amount"})
    df_list.append(df)
del file, df

df_dict = dict(zip(spentmed_file_names, df_list)) #excellent. Now we have all our files and their file names

#Import metadata; will need later
metadata = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/Metadata_for_human_and_microbial_communities_Burns_data_V2.csv")
metadata_slimmed = metadata[["Patient_Blind_ID", "host_sampleid", "host_description"]]
metadata_slimmed.drop_duplicates(subset=["host_sampleid"], inplace=True) #get a copy of slice error. Whatever
metadata_slimmed["Patient_Blind_ID"]=metadata_slimmed["Patient_Blind_ID"].apply(str) #set PBID to string; needed for merging later
metadata_slimmed["match_col"] = metadata_slimmed["Patient_Blind_ID"]+"_"+metadata_slimmed["host_description"] #create match col for merging

#now, to extract (and add back in) metadata from the filenames
edited_df_list = []
for key, value in df_dict.items():
    working_df = value
    orig_string = key
    filename_list = orig_string.split("_")
    pbid_string = filename_list[1] #pick out the element that represents the PIBD
    description_string = filename_list[2] #pick out the element that represents the description (tumor/normal)
    working_df["Patient_Blind_ID"] = str(pbid_string)
    working_df["Description"] = description_string
    working_df["match_col"] = filename_list[1]+"_"+filename_list[2] 
    working_df2= pd.merge(left=working_df, right= metadata_slimmed, how='left', left_on="match_col", right_on="match_col")
    edited_df_list.append(working_df2)
del working_df, orig_string, filename_list, pbid_string, description_string, working_df2    

#now, stack all DFs vertically into a single DF
concat_spentmed = pd.concat(edited_df_list, axis=0) #concatenate
concat_spentmed = concat_spentmed.drop(columns={"Patient_Blind_ID_y", "Description", "match_col"})

#GOOD! now, save
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_analyses")
concat_spentmed.to_csv("04.01.22_combining_host_and_microbiomes_EUstd_host_spent_medium.csv")

#%reset clear out variables

#############PART 2: MICOM cell growth rates
#let's start by bringing in some files
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_results/EUstd_CORDA_to_MICOM_part2_results/")

comm_GR_file_list = glob.glob('*_EUstd_GRs.csv')
comm_GR_file_names = []

for elem in comm_GR_file_list:
    sample_name = re.sub("_EUstd_GRs.csv", "", elem)  
    comm_GR_file_names.append(sample_name)        
del elem, sample_name    

comm_GRs =  {}
for elem in comm_GR_file_list:
    output_file = pd.read_csv(elem)
    sample_name = re.sub("_EUstd_GRs.csv", "", elem) #delete extra crap on string.
    for column in output_file:
        costring1 = re.sub("<CommunitySolution ", "", column) #delete first part of string from solution
        costring2 = costring1[:-19] #delete last 19 characters of string. These all differ, so can't match text/symbols
        comm_GRs.update({sample_name: costring2})
del elem, column, costring1, costring2, output_file, sample_name

#this is all previously borrowed. Now, the trick is to extract the metadata from the sample name.
comm_GRs_df = pd.DataFrame(comm_GRs.items(), columns=['file_name', 'growth_rate']) #first, convert comm_GRs to a dataframe
comm_GRs_df = comm_GRs_df.set_index("growth_rate") #reset index (needed because growth rate will otherwise be lost)
comm_GRs_df = comm_GRs_df['file_name'].str.split('_', expand=True) #then, break the file_name column into multiple columns based on delimiter (here, _)
comm_GRs_df.reset_index(inplace=True) #reset index to return growth_rate to a column value
#need to merge columns 7 and 8 for microbiome sample id
comm_GRs_df["microbiome_id"]=  comm_GRs_df[2] +'_'+ comm_GRs_df[3]
#Now, recombine columns into useful headings
comm_GRs_df.rename(columns={1: 'Patient_Blind_ID', 4: "microbiome_description", 7: "host_id", 8: "host_description"}, inplace= True)
comm_GRs_df.drop([0,2,3,5,6,9], axis=1, inplace=True) #drop extra columns
#let's also add a new column, combo_type, to help us with our graphing
comm_GRs_df["combo_type"] = comm_GRs_df["host_description"]+"_"+comm_GRs_df["microbiome_description"]
#GOOD! now, save
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_analyses")
comm_GRs_df.to_csv("04.01.22_combining_host_and_microbiomes_EUstd_microbiome_comm_GRs.csv")

del comm_GR_file_list, comm_GR_file_names, comm_GRs, comm_GRs_df

#Part 3: MICOM indiv spp growth rates
#let's start by bringing in some files
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_results/EUstd_CORDA_to_MICOM_part2_results/")

indiv_GR_file_list = glob.glob('*_EUstd_indiv_spp_GRs.csv')

df_list = []
for file in indiv_GR_file_list:
    df= pd.read_csv(file)
    df=df.rename(columns={"compartment": "OTU_ID"})
    df["combo_type"] = df["host_model_description"]+"_"+df["microbiome_model_description"]
    df_list.append(df)
del file, df

concat_indiv_spp_GRs = pd.concat(df_list, axis=0) #concatenate
#dang, having those metadata columns really makes life easier.

#GOOD! now, save
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_analyses")
concat_indiv_spp_GRs.to_csv("04.01.22_combining_host_and_microbiomes_EUstd_microbiome_indiv_spp_GRs.csv")

del concat_indiv_spp_GRs, df_list, indiv_GR_file_list

##############Part 4: MICOM fluxes
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_results/EUstd_CORDA_to_MICOM_part2_results/")
flux_file_list = glob.glob('*_EUstd_fluxes.csv')

sample_names_list = []
for file in flux_file_list:
    sample_name = re.sub("_EUstd_fluxes.csv", "", file) #delete extra crap on string.  
    sample_names_list.append(sample_name)
del file, sample_name

df_list = []
for file in flux_file_list:
    df= pd.read_csv(file)
    df=df.rename(columns={"compartment": "OTU_ID"})
    df_list.append(df)
del file, df

df_dict = dict(zip(sample_names_list, df_list)) #add back later, for transposed DFs

######step 2: clean up files and prepare for merging
df_list_cleaned = []
for key, value in df_dict.items():
    df_cleaned = value.set_index('OTU_ID')
    df_cleaned= df_cleaned.drop("medium", axis=0)
    df_transposed = df_cleaned.T
    df_transposed.reset_index(inplace=True)
    df_transposed = df_transposed.rename(columns={"index": "fluxpath_id"})
    df_long =pd.melt(df_transposed,id_vars=['fluxpath_id'],var_name='otu_id', value_name='fluxpath_values')
    df_long['SampleID']=str(key)
    tempsplit = df_long['SampleID'].str.split('_', expand=True) #then, break the file_name column into multiple columns based on delimiter (here, _)
    df_long["Patient_Blind_ID"] = tempsplit[1]
    df_long["microbiome_id"] = tempsplit[2] +'_'+ tempsplit[3]
    df_long["microbiome_description"] = tempsplit[4]
    df_long["host_id"] = tempsplit[7]
    df_long["host_description"] = tempsplit[8]
    df_list_cleaned.append(df_long)
del key, value, df_cleaned, df_transposed, df_long    


concat_flux = pd.concat(df_list_cleaned, axis=0) #concatenate

#GOOD! now, save
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_analyses")
concat_flux.to_csv("04.01.22_combining_host_and_microbiomes_EUstd_microbiome_fluxes.csv")

del concat_flux, df_list, flux_file_list













































