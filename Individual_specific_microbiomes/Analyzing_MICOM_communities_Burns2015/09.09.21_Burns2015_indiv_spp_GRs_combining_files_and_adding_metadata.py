#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 13:51:22 2021

@author: adamo010
"""

import scipy
import numpy as np
import csv
import subprocess
import os
import pandas as pd
import fnmatch
import shutil
import glob
import re
import copy
from functools import reduce

#the purpose fo this code is to create a single csv file containing individual species growth rates from all tradeoff parameters run. 
#Also, add taxonomic information because that's interesting (or at least more so than otu_id).
#based on 02.04.21_combining_species_specific_GRs_over_TD_values.py and 
#Date_TBA_combining_species_specific_GRs_from_MSI_runs.py and
#08.31.21_Hale2018_indiv_spp_GRs_combining_files_and_adding_metadata.py

##################step 1: import files
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Common_growth_medium_runs/June2021_rerunning_MICOM_on_MSI_with_fluxes_on/Sept2021_V30_output")

#start by importing all the files. 
filelist = glob.glob('*_V30_Burns_fluxes_on_indiv_spp_GRs.csv')

#make a list of sample names
sample_names_list = []
for file in filelist:
    sample_name = re.sub("_V30_Burns_fluxes_on_indiv_spp_GRs.csv", "", file) #delete extra crap on string.
    sample_names_list.append(sample_name)
del file, sample_name    

#make a list of dataframes
df_list = []
for file in filelist:
    df= pd.read_csv(file)
    df_list.append(df)
del file, df

#zip the filename list and dataframe list into a dictionary. 
df_dict = dict(zip(sample_names_list, df_list))

##################step 2: combine all samples into single dataframe and add taxonomy data
#import taxonomy information as preprocessed_taxonomy
#note that this is the FULL taxonomy table (i.e. contaminants not removed.)
preprocessed_taxonomy = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Burns_2015_OTU_id_taxonomy_table_tax_levels_separated.csv")
preprocessed_taxonomy.drop("Unnamed: 0", axis=1, inplace=True)

dataframe_list = {} #create a dictionary of growth rate files merged with taxonomies

def combining_across_samples(sample_id):
    for file in os.listdir():
        if sample_id in str(file) and str("V30_Burns_fluxes_on_indiv_spp_GRs") in str(file):
            input_file = pd.read_csv(file)
            input_file["growth_rate"][input_file["growth_rate"] < 0.000001] = 0 #sets all values less than 10-6 as 0
            input_plus_taxonomy = pd.merge(left= input_file, right=preprocessed_taxonomy, how="left", left_on="compartments", right_on="OTU_ID" )
            input_plus_taxonomy.drop("compartments", axis=1, inplace=True)
            input_plus_taxonomy["Sample_ID"] = str(sample_id)
            dataframe_list[str(sample_id)] = input_plus_taxonomy
    return
        
for sample in sample_names_list:
    combining_across_samples(sample)
del sample    

combined_indiv_GRs = pd.concat(dataframe_list.values(), ignore_index=True) #this has not dropped taxa with 0 growth rates.

##################step 3: import and add metadata
metadata = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Burns_2015_sample_metadata.csv")#wtf, why are therea bunch of extra rows here?
metadata.dropna(subset = ["Patient_Blind_ID"], inplace=True)

combined_indiv_GRs_with_metadata = pd.merge(left=combined_indiv_GRs, right=metadata, how='left', left_on='Sample_ID', right_on='SampleID')
combined_indiv_GRs_with_metadata= combined_indiv_GRs_with_metadata.drop(columns=['SampleID'])

##################step 4: collapsing to get a difference in growth rates between tumor and normal samples.
#all right, here are a couple of helpful lists for sorting:
tumor_samples = []
normal_samples = []
temp_df= metadata.loc[metadata['Description'] == 'tumor']
for row in temp_df["SampleID"]:
    tumor_samples.append(row)
temp_df= metadata.loc[metadata['Description'] == 'normal']
for row in temp_df["SampleID"]:
    normal_samples.append(row)
del temp_df, row  

#create a column with a unique value for every otu_id/patient_blind_id combo. Can use these to merge on when we reunite these.
combined_indiv_GRs_with_metadata["unique_id"]= combined_indiv_GRs_with_metadata["OTU_ID"].astype(str) + combined_indiv_GRs_with_metadata["Patient_Blind_ID"].astype(str)

#now, create separate dataframes for each site (tumor/normal)
combined_indiv_GRs_with_metadata_tumor =  combined_indiv_GRs_with_metadata[combined_indiv_GRs_with_metadata['Description'] == "tumor"] 
combined_indiv_GRs_with_metadata_normal =  combined_indiv_GRs_with_metadata[combined_indiv_GRs_with_metadata['Description'] == "normal"] 

#now, rename each growth_rate column as growth_rate_normal or growth_rate_tumor
combined_indiv_GRs_with_metadata_tumor = combined_indiv_GRs_with_metadata_tumor.rename(columns={'growth_rate': 'growth_rate_tumor'})
combined_indiv_GRs_with_metadata_normal = combined_indiv_GRs_with_metadata_normal.rename(columns={'growth_rate': 'growth_rate_normal'})

#I think I need to strip all the columns off the tumor (smaller) OTUs because otherwise everything is a mess
#tumor is still smaller in Hale dataset, which is reassuring
combined_indiv_GRs_with_metadata_tumor = combined_indiv_GRs_with_metadata_tumor[['growth_rate_tumor','unique_id']]

#okay NOW try this join thing.
combined_indiv_GRs_left = pd.merge(left=combined_indiv_GRs_with_metadata_normal, right=combined_indiv_GRs_with_metadata_tumor, how='left', left_on='unique_id', right_on='unique_id')
combined_indiv_GRs_inner = pd.merge(left=combined_indiv_GRs_with_metadata_normal, right=combined_indiv_GRs_with_metadata_tumor, how='inner', left_on='unique_id', right_on='unique_id')

#so... Now I will subtract these columns. I imagine it won't work with the left joins, because you can't subtract NANs, so I'll just do it with the inner ones
combined_indiv_GRs_left["growth_rate_difference"] = combined_indiv_GRs_left["growth_rate_tumor"] - combined_indiv_GRs_left["growth_rate_normal"]
combined_indiv_GRs_inner["growth_rate_difference"] = combined_indiv_GRs_inner["growth_rate_tumor"] - combined_indiv_GRs_inner["growth_rate_normal"]

#just out of curiosity, if I drop the nans in point_X_left, does it make point_X_inner?
demo = combined_indiv_GRs_left.copy(deep=True)
demo.dropna(subset=['growth_rate_difference'], inplace=True)
#yeah, dog. Same number of rows as combined_indiv_GRs_inner.
del demo

##################step 5: collapsing by genus and family. This output is used for graphing boxplots!
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/") #analysis directory
#collapsing by genus
temp = copy.deepcopy(combined_indiv_GRs)
temp.drop(["species"], axis=1, inplace=True) #not dropping taxonomy_orig, b/c we don't have it. Capitalization also an issue here
collapsed_GRs = temp.groupby(["OTU_ID", "Sample_ID", 'kingdom', 'phylum', 'class', 'order', 'family', 'genus']).mean()
collapsed_GRs.columns = ['mean_genus_GR']
collapsed_GRs = collapsed_GRs.reset_index()
collapsed_GRs_with_metadata = pd.merge(left=collapsed_GRs, right=metadata, how='left', left_on='Sample_ID', right_on='SampleID')
collapsed_GRs_with_metadata.drop(["Sample_ID"], axis=1, inplace=True)
collapsed_GRs_with_metadata.to_csv("09.09.21_Burns2015_indiv_spp_GRs_collapsed_to_genus.csv")
del temp

#collapsing by family
temp2 = copy.deepcopy(combined_indiv_GRs)
temp2.drop(["species", "genus"], axis=1, inplace=True)
collapsed_GRs2 = temp2.groupby(["OTU_ID", "Sample_ID", 'kingdom', 'phylum', 'class', 'order', 'family']).mean()
collapsed_GRs2.columns = ['mean_family_GR']
collapsed_GRs2 = collapsed_GRs2.reset_index()
collapsed_GRs_with_metadata2 = pd.merge(left=collapsed_GRs2, right=metadata, how='left', left_on='Sample_ID', right_on='SampleID')
collapsed_GRs_with_metadata2.drop(["Sample_ID"], axis=1, inplace=True)
collapsed_GRs_with_metadata2.to_csv("09.09.21_Burns2015_indiv_spp_GRs_collapsed_to_family.csv")
del temp2

##################step 6: collapsed differences: use combined_indiv_GRs_inner
temp = copy.deepcopy(combined_indiv_GRs_inner)
temp.drop(["taxonomy_orig", "species", "Description", "Patient_Blind_ID", "Site", "MSI_status", "Stage",
           "unique_id"], axis=1, inplace=True)
collapsed_GRs_diffs = temp.groupby(["OTU_ID", "Sample_ID", 'kingdom', 'phylum', 'class', 'order', 'family', 'genus']).mean()
#add an extra step to re-calculate growth rate difference; shouldn't be averaging differences. 
collapsed_GRs_diffs.columns = ['mean_genus_GR_normal', "mean_genus_GR_tumor", "mean_GR_difference"]
#collapsed_GRs_point_one["mean_GR_differences"] = collapsed_GRs_point_one["mean_genus_GR_tumor"] - collapsed_GRs_point_one["mean_genus_GR_normal"]
collapsed_GRs_diffs = collapsed_GRs_diffs.reset_index()
collapsed_GRs_diffs_with_metadata = pd.merge(left=collapsed_GRs_diffs, right=metadata, how='left', left_on='Sample_ID', right_on='SampleID')
collapsed_GRs_diffs_with_metadata.drop(["Sample_ID"], axis=1, inplace=True)
collapsed_GRs_diffs_with_metadata.to_csv("09.09.21_Burns2015_indiv_spp_GRs_collapsed_to_genus_diffs.csv")
