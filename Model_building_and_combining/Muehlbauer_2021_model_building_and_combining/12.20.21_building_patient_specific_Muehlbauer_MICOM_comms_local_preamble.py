#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 09:28:57 2021

@author: adamo010
"""

#somewhat based on 08.10.20_building_patient_specific_MICOM_comms.py
#very based on 08.05.21_building_patient_specific_Hale_MICOM_comms_localtest.py
#and 10.22.21_building_patient_specific_Niccolai_MICOM_comms_local_preamble.py

#The purpose of this code is to generate .pickle MICOM community files for each patient in
#the Niccolai dataset (80 communities total, 2 sites per 40 patients).

#The general outline/plan is as follows:
#1. Import the abundances and the otu name/model name file
#2. Clip taxonomy off abundances (front and back of table)
#3. Match model file names to each otu
#4. Create new dataframe where only a given (sample) column is included and only non-zero abundance rows are included.
#5. Save that dataframe as the Sample_number
#6. Import all those and make the communities
#7. Save the communities as .pickle files named by Sample_number\

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
import pickle
from ast import literal_eval

#I realize that I'll need to re-do abundance filtering based on the OTU table that I actually ended up with
#bring in absolute abundance OTU table and do relative abundances from a filtered version of that.
#importing files

os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Muehlbauer_2021_model_picking/")

#import info for all models; these were geneated in 12.03.21_combining_all_models_for_Muehlbauer2021_OTUs.py
agora_model_otu_data = pd.read_pickle('Muehlbauer2021_AGORA_models_that_need_OTU_builds.pickle')
carveme_model_otu_data = pd.read_pickle('Muehlbauer2021_CarveMe_models_that_need_otu_builds.pickle')
mambo_model_otu_data = pd.read_pickle('Muehlbauer2021_MAMBO_models_that_need_OTU_builds.pickle')

#build a new dataframe that looks like 06.22.20_summary_of_otus_with_models.csv
#start with each database-specific otu, then merge them
agora_slimmed = agora_model_otu_data[["OTU_ID", "all_agora_model_file_names"]].copy()
mambo_slimmed = mambo_model_otu_data[["OTU_ID", "all_mambo_matched_models"]].copy()
carveme_slimmed = carveme_model_otu_data[["OTU_ID", "all_carveme_matched_models"]].copy()

#add a column that refers to the database from which the OTU's models came
agora_slimmed['model_db']='agora'
mambo_slimmed['model_db'] = 'mambo'
carveme_slimmed['model_db']='carveme'

#count the number of models in each row
def counting_models(df, col_name):
    counting_col = []
    for row in df[col_name]:
        num_models = len(row)
        counting_col.append(num_models)
    df["number_of_models"] = counting_col
    return
counting_models(agora_slimmed, "all_agora_model_file_names")
counting_models(mambo_slimmed, "all_mambo_matched_models")
counting_models(carveme_slimmed, "all_carveme_matched_models")

#add a new column called final_model_name: 
def otu_file_names(df, col_name):
    file_names_list = []
    for row in df[col_name]:
        file_name = str(row)+"_model.xml"
        file_names_list.append(file_name)
    df["final_model_name"]= file_names_list
    return
otu_file_names(agora_slimmed, "OTU_ID")    
otu_file_names(mambo_slimmed, "OTU_ID")    
otu_file_names(carveme_slimmed, "OTU_ID")    
        
#rename all the model_file_names colums so they are consistent. new-old
agora_slimmed.rename(columns = {'all_agora_model_file_names': 'model_file_names'}, inplace=True)
mambo_slimmed.rename(columns = {'all_mambo_matched_models': 'model_file_names'}, inplace=True)
carveme_slimmed.rename(columns = {'all_carveme_matched_models': 'model_file_names'}, inplace=True)

#now, the great merge
all_otus_with_models = pd.concat([agora_slimmed, mambo_slimmed, carveme_slimmed])
#SLAM. DUNK. 1489 OTUs. 
all_otus_with_models.to_csv("12.20.21_summary_of_Muehlbauer_OTUs_with_models.csv")

#cool. now need an OTU table that corresponds to all these OTUs. 
#don't actually have one, so do some fenablgilng.
otus_with_low_tax_ids=pd.read_csv("12.01.21_Muehlbauer2021_OTUs_with_low_tax_ids.csv")
prefiltered_otus = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal/12.01.21_Muehlbauer2021_asv_table_abundance_filtered_0.1_any_sample_cutoff.csv") 
#make a list of low tax OTUs to remove from prefiltered otus
low_tax_otus = otus_with_low_tax_ids["OTU_ID"].tolist()
#filter prefiltered_otus to remove low_tax_ids
model_matched_otus = prefiltered_otus[~prefiltered_otus["OTU_ID"].isin(low_tax_otus)]
#NEW step: filter out OTUs for which we could not build models (there are five, which I shall simply list here, as they have not been
#pulled out and stuck in a file.)
unmatched_otus = ["ASV_134", "ASV_16", "ASV_71", "ASV_334", "ASV_326"] 
model_matched_otus = model_matched_otus[~model_matched_otus["OTU_ID"].isin(unmatched_otus)]

#clean up
del otus_with_low_tax_ids, prefiltered_otus, agora_model_otu_data, agora_slimmed, carveme_model_otu_data, carveme_slimmed, low_tax_otus
del mambo_slimmed, mambo_model_otu_data, unmatched_otus

###############OKAY, preamble stuff done. 
#first, generate a list of sample names. each of these needs a micom community generated.
sample_names = []            #this is a list of all the column names that will need to be scanned for values >0.001
for col in model_matched_otus.columns:       
    sample_names.append(str(col))
del col
sample_names.remove("taxonomy_orig")
sample_names.remove("OTU_ID")
sample_names.remove('keep_or_toss')
#okay, now we have a list 12 sample names.

#next step is to add a column to model_matched_otus that lists that otu's model_id
model_matched_otus["OTU_ID"] = model_matched_otus["OTU_ID"].astype(str) #convert otu_ids to strings
all_otus_with_models["OTU_ID"] = all_otus_with_models["OTU_ID"].astype(str)
model_file_dict = dict(zip(all_otus_with_models.OTU_ID, all_otus_with_models.final_model_name))
model_matched_otus['otu_file_name']= model_matched_otus["OTU_ID"].apply(lambda x: model_file_dict.get(x))

#the next steps I want to do are:
#1) create new dataframes, one for each sample
#2) for each dataframe, only keep non-zero abundance rows
#3) for each dataframe, only keep the following columns: taxonomy, OTU_id, model_file_name, Sample_XX(relative abundance))

#probably have to create a function for this.
#I KNOW you're not supposed to iteratively create variables, but I'm doing it anyway. 

#ALSO NOTE: probably have to create a new directory before running the below function

def subsetting_OTU_table(sample_name):
    subset_df = model_matched_otus[model_matched_otus[str(sample_name)] > 0]
    #the above line preserves all rows (OTUs) where GR >0 in the column corresponding to sample_name in OTU_matched_models
    subset_df2 = subset_df[["taxonomy_orig", "OTU_ID", "otu_file_name", str(sample_name)]] #drop extra columns.
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/Muehlbauer_2021_model_building_and_combining/Muehlbauer_OTU_specific_abundance_tables/")
    subset_df2.to_csv("_abundances_and_models.csv")
    for file in os.listdir():                   
        src=file
        if fnmatch.fnmatch(file, "_abundances_and_models.csv"):
            dst = str(sample_name)+file
            os.rename(src,dst)
    return

subsetting_OTU_table('SRR12572550') #this was a test

#now we do the whole one:
for elem in sample_names:
    subsetting_OTU_table(str(elem))
del elem    

#fucking incredible, dude. 
    
#we now have a set of tables containing the necessary information for MICOM. These files are all named as follows:
#sampleID_abundances_and_models.csv
#Now we will import these and make communities with them. Hopefully on MSI. 

#but first, let's make a text file that is a list of these.
a= open("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/Muehlbauer_2021_model_building_and_combining/Muehlbauer2021_OTU_abundances_all.txt", "a")
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/Muehlbauer_2021_model_building_and_combining/Muehlbauer_OTU_specific_abundance_tables/")
for file in os.listdir():
    filename = str(file)
    a.write(str(filename) + os.linesep)
    
#big W for this Friday.

#prepare a list of pickle files for running MICOM on MSI    
sample_names_list = []
for item in sample_names:
    pickle_name = str(item)+"_community.pickle"
    sample_names_list.append(pickle_name)
    
os.chdir('/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Running_MICOM_communities_Muehlbauer2021/')

with open('Muehlbauer2021_pickle_file_list.txt', 'w') as f:
    for item in sample_names_list:
        f.write("%s\n" % item)













