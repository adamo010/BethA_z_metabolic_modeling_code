#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 11:31:21 2021

@author: adamo010
"""
#somewhat based on 08.10.20_building_patient_specific_MICOM_comms.py

#The purpose of this code is to generate .pickle MICOM community files for each patient in
#the Hale dataset (88 communities total, 2 sites per 44 patients).

#The general outline/plan is as follows:
#1. Import the abundances and the otu name/model name file
#2. Clip taxonomy off abundances (front and back of table)
#3. Match model file names to each otu
#4. Create new dataframe where only a given (sample) column is included and only non-zero abundance rows are included.
#5. Save that dataframe as the Sample_number
#6. Import all those and make the communities
#7. Save the communities as .pickle files named by Sample_number\

#note that this takes place in the MICOM_CRC_FBA_02.2020 environment
#the directory I'm going to set as home is /Users/adamo010/Documents/MICOM_CRC_FBA
#the output files shoudl go into the following directory: /Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes

#files that I will need: 
#06.22.20_summary_of_otus_with_models.csv.
#06.16.20_otu_table_abundance_filtered_0.1.csv

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

#The August2020 version of this script used a file called 06.22.20_summary_of_otus_with_models.csv
#I don't have the equivalent of that, so make it myself here.
#use code from 06.22.20_summarizing_model_assignment_results.py and model-otu info made 
#in 07.29.21_combining_all_models_for_Hale2018_otus.py

os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Hale_2018_model_picking/")

#import info for all models
agora_model_otu_data = pd.read_pickle('AGORA_models_that_need_OTU_builds.pickle')
mambo_model_otu_data = pd.read_pickle('MAMBO_models_that_need_OTU_builds.pickle')
carveme_model_otu_data = pd.read_pickle('CarveMe_models_that_need_OTU_builds.pickle')

#build a new dataframe that looks like 06.22.20_summary_of_otus_with_models.csv
#start with each database-specific otu, then merge them
agora_slimmed = agora_model_otu_data[["OTU_ID", "OTU_proxy_ID", "all_agora_model_file_names"]].copy()
mambo_slimmed = mambo_model_otu_data[["OTU_ID", "OTU_proxy_ID", "all_mambo_matched_models_x"]].copy()
carveme_slimmed = carveme_model_otu_data[["OTU_ID", "OTU_proxy_ID", "all_carveme_matched_models"]].copy()

#add a column that refers to the database from which the OTU's models came
agora_slimmed['model_db']='agora'
mambo_slimmed['model_db']='mambo'
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
counting_models(mambo_slimmed, "all_mambo_matched_models_x")
counting_models(carveme_slimmed, "all_carveme_matched_models")

#add a new column called final_model_name: 
def otu_file_names(df, col_name):
    file_names_list = []
    for row in df[col_name]:
        file_name = str(row)+"_model.xml"
        file_names_list.append(file_name)
    df["final_model_name"]= file_names_list
    return
otu_file_names(agora_slimmed, "OTU_proxy_ID")    
otu_file_names(mambo_slimmed, "OTU_proxy_ID")    
otu_file_names(carveme_slimmed, "OTU_proxy_ID")    
        
#rename all the model_file_names colums so they are consistent. new-old
agora_slimmed.rename(columns = {'all_agora_model_file_names': 'model_file_names'}, inplace=True)
mambo_slimmed.rename(columns = {'all_mambo_matched_models_x': 'model_file_names'}, inplace=True)
carveme_slimmed.rename(columns = {'all_carveme_matched_models': 'model_file_names'}, inplace=True)

#now, the great merge
all_otus_with_models = pd.concat([agora_slimmed, mambo_slimmed, carveme_slimmed])
#SLAM. DUNK. 1489 OTUs. 
all_otus_with_models.to_csv("08.05.21_summary_of_Hale_OTUs_with_models.csv")

#cool. now need an OTU table that corresponds to all these OTUs. 
#don't actually have one, so do some fenablgilng.
otus_with_low_tax_ids=pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Hale_2018_model_picking/07.13.21_Hale2018_OTUs_with_low_tax_ids.csv")
prefiltered_otus = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal/07.13.21_Hale2018_otu_table_abundance_filtered_0.1_any_sample_cutoff.csv") 
#make a list of low tax OTUs to remove from prefiltered otus
low_tax_otus = otus_with_low_tax_ids["OTU_ID"].tolist()
#filter prefiltered_otus to remove low_tax_ids
model_matched_otus = prefiltered_otus[~prefiltered_otus["OTU_ID"].isin(low_tax_otus)]
#clean up
del otus_with_low_tax_ids, prefiltered_otus, agora_model_otu_data, agora_slimmed, carveme_model_otu_data, carveme_slimmed,\
    mambo_model_otu_data, mambo_slimmed, low_tax_otus
    
###############OKAY, preamble stuff done. 
#first, generate a list of sample names. each of these needs a micom community generated.
sample_names = []            #this is a list of all the column names that will need to be scanned for values >0.001
for col in model_matched_otus.columns:       
    sample_names.append(str(col))
del col
sample_names.remove("taxonomy_orig")
sample_names.remove("OTU_ID")
sample_names.remove('keep_or_toss')
#okay, now we have a list of 92 sample names.

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

def subsetting_OTU_table(sample_name):
    subset_df = model_matched_otus[model_matched_otus[str(sample_name)] > 0]
    #the above line preserves all rows (OTUs) where GR >0 in the column corresponding to sample_name in OTU_matched_models
    subset_df2 = subset_df[["taxonomy_orig", "OTU_ID", "otu_file_name", str(sample_name)]] #drop extra columns.
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Building_MICOM_communities_Hale2018/Hale_OTU_specific_abundance_tables")
    subset_df2.to_csv("_abundances_and_models.csv")
    for file in os.listdir():                   
        src=file
        if fnmatch.fnmatch(file, "_abundances_and_models.csv"):
            dst = str(sample_name)+file
            os.rename(src,dst)
    return

#subsetting_OTU_table('s006F13xB1H04') #this was a test

#now we do the whole one:
for elem in sample_names:
    subsetting_OTU_table(str(elem))
del elem    
    
#we now have a set of tables containing the necessary information for MICOM. These files are all named as follows:
#sampleID_abundances_and_models.csv
#Now we will import these and make communities with them (and very probably break my poor computer).

#let's start with Sample 27.
s006F13xB1H04_df = pd.read_csv("s006F13xB1H04_abundances_and_models.csv")    
s006F13xB1H04_df.drop(["Unnamed: 0", "taxonomy_orig"], axis=1, inplace=True)    #get rid of extra column.
#now why tf did I keep the taxonomy column if we're not even going to use it?
#anyway.
       
from micom import Community    
#sample_27_community = Community(Sample_27_df, solver='gurobi')    
#can't do this raw; need an id column
#christ, what a horrible way to phrase it. What I"m doing here is renaming columns to match those in the micom Community function
s006F13xB1H04_df.rename(columns={'OTU_ID': 'id', 's006F13xB1H04': 'abundance', 'otu_file_name': 'file'}, inplace=True)
#I keep getting this "AttributeError: Can only use .str accessor with string values!" error. Maybe the stupid
#program is reading my id column as not strings?
#df['DataFrame Column'] = df['DataFrame Column'].astype(float)
s006F13xB1H04_df['id']= s006F13xB1H04_df["id"].astype(str)
#so, that worked, but now need to change directory to /Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/OTU_models
#there is probably a way to do this in the command, but I don't know it. 
s006F13xB1H04_community = Community(s006F13xB1H04_df, solver='gurobi')     
#and now we wait. Looking at 10-20 seconds per model, this particular sample has 58 OTUs, so... probably about 20 minutes.
#yep. About 22 minutes.
#now, save to pickle.
s006F13xB1H04_community.to_pickle("s006F13xB1H04_community.pickle")
#cool.

#in the meantime, I'll work on writing parallelization code.

#AT this point, we can move onto MSI. Write the code here, but copy the below stuff to an MSI file

################interm: running MAMBO locally
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Hale_2018_model_picking/")
all_mambo_merged_condensed2 = pd.read_pickle("MAMBO_models_that_need_OTU_builds.pickle")
all_mambo_merged_condensed2= all_mambo_merged_condensed2.drop(columns=["OTU_ID"])
#interm- pull out AGORA models.
#indices 7, 28, 62, 63 OR where .xml is present
weird_extra_agora_otus= all_mambo_merged_condensed2.iloc[[7,28,62,63]]
all_mambo_merged_condensed3 = all_mambo_merged_condensed2.drop([7,28,62,63])

mambo_model_file_dict = dict(zip(all_mambo_merged_condensed3.OTU_proxy_ID,all_mambo_merged_condensed3.all_mambo_matched_models_x.apply(lambda x: literal_eval(str(x)))))    
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/MAMBO_models/")
counter = 0
for key, value in mambo_model_file_dict.items():
    filename = str(key)+"_model.xml"
    otu_model = micom.util.join_models(value, str(key))
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/Hale_2018_model_building_and_combining/Hale_2018_otu_specific_models/")
    cobra.io.write_sbml_model(otu_model, filename)
    del otu_model
    counter= counter+1
    print(str(counter))
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/MAMBO_models/")
del counter, filename, key, value
#rerun agora extra crap too.
####import agora model key
agora_model_file_key= pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/agora_model_id_and_file_key.csv",\
                                  index_col=0, squeeze=True).to_dict()
matched_model_id_list = weird_extra_agora_otus["all_mambo_matched_models_x"].to_list() #convert model ID column to list
all_model_files = [] 
#for each model id in the list of model_ids associated with each row (OTU) in the dataframe, match with a model file
#from the agora_model_file_key dictionary
for elem in matched_model_id_list:
    model_file_list = []
    for subelem in elem:
        model_file_string = subelem.replace("_"," ")
        #print(model_file_string)
        matched_model_file_string = str(agora_model_file_key[model_file_string])
        #print(matched_model_file_string)
        model_file_list.append(matched_model_file_string)
    all_model_files.append(model_file_list)
del elem, subelem, model_file_string  
#now, have a list called all_model_files that can be stuck on as a column in the agora dataframe
weird_extra_agora_otus['all_agora_model_file_names_fixed'] = all_model_files
#now, make AGORA models here also.
weird_agora_model_file_dict = dict(zip(weird_extra_agora_otus.OTU_proxy_ID,weird_extra_agora_otus.all_agora_model_file_names_fixed.apply(lambda x: literal_eval(str(x)))))    
counter = 0
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/AGORA_models/")
for key, value in weird_agora_model_file_dict.items():
    filename = str(key)+"_model.xml"
    otu_model = micom.util.join_models(value, str(key))
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/Hale_2018_model_building_and_combining/Hale_2018_otu_specific_models/")
    cobra.io.write_sbml_model(otu_model, filename)
    del otu_model
    counter= counter+1
    print(str(counter))
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/AGORA_models/")
del counter, filename, key, value  
    
###############interm over

#first, create a list of all the files with each sample's OTU table    
sample_OTU_table_files = []
for file in glob.glob("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Building_MICOM_communities_Hale2018/Hale_OTU_specific_abundance_tables/*abundances_and_models.csv"):
    sample_OTU_table_files.append(file)   
del file
sample_OTU_table_samplenames = []
for elem in sample_OTU_table_files:
    elem = elem.replace("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Building_MICOM_communities_Hale2018/Hale_OTU_specific_abundance_tables/","")
    elem = elem.replace("_abundances_and_models.csv","")
    sample_OTU_table_samplenames.append(elem)
del elem

#merge these into a dictonary; then we can iterate over both sets.
sample_dict = {sample_OTU_table_samplenames[i]: sample_OTU_table_files[i] for i in range(len(sample_OTU_table_samplenames))}

def creating_MICOM_tables(sample_file):
    file_name = str(sample_file)
    sample_name = file_name.replace("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Building_MICOM_communities_Hale2018/Hale_OTU_specific_abundance_tables/","")
    sample_name = sample_name.replace("_abundances_and_models.csv","")
    sample_name.df = pd.read_csv("Individual_specific_microbiomes/sample_file")
    return

#step 1: importing and editing all dataframes
edited_sample_dict = {}
for key, value in sample_dict.items():
    dataframe = pd.read_csv(value)
    dataframe.rename(columns = {'OTU_ID': 'id', key: 'abundance', 'otu_file_name': 'file'}, inplace=True)
    dataframe['id']= dataframe['id'].astype(str)
    edited_sample_dict.update( {key : dataframe})
del key, value    

#this is not ideal, but... yeah. all values are called dataframe because APPARENTLY dynamically creating variables
#is wrong.

#step 2: running MICOM
#oh god, this is going to break my machine.
finished_communities = {}
for key, value in edited_sample_dict.items():
    sample_community = Community(value, solver='gurobi')
    finished_communities.update( {key: sample_community})     
#this took, I shit you not, five days. 

#step 3: saving as .pickle
for key, value in finished_communities.items():
  value.to_pickle("_community.pickle")
  for file in os.listdir():
    src= file
    if fnmatch.fnmatch(file, "_community.pickle"):
        dst = str(key)+file
        os.rename(src, dst)    
#phew, that took a few days too, didn't it? Let's see if they all worked.


###################fixing MSI runs for AGORA models.
#%reset #type this into the console to clear variables.
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Hale_2018_model_picking/")
all_agora_merged_condensed2 = pd.read_pickle("AGORA_models_that_need_OTU_builds.pickle")
all_agora_merged_condensed2= all_agora_merged_condensed2.drop(columns=["OTU_ID"])

#now, generate models  
agora_model_file_dict = dict(zip(all_agora_merged_condensed2.OTU_proxy_ID,all_agora_merged_condensed2.all_agora_model_file_names.apply(lambda x: literal_eval(str(x)))))    
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/AGORA_models/")
#break agora_model_file_dict into parts
agora_dict_A = dict(list(agora_model_file_dict.items())[len(agora_model_file_dict)//2:])
agora_dict_B = dict(list(agora_model_file_dict.items())[:len(agora_model_file_dict)//2])
agora_dict_1 = dict(list(agora_dict_A.items())[len(agora_dict_A)//2:])
agora_dict_2 = dict(list(agora_dict_A.items())[:len(agora_dict_A)//2])
agora_dict_3 = dict(list(agora_dict_B.items())[len(agora_dict_B)//2:])
agora_dict_4 = dict(list(agora_dict_B.items())[:len(agora_dict_B)//2])
del agora_dict_A, agora_dict_B

counter = 0
for key, value in agora_dict_1.items():
    filename = str(key)+"_model.xml"
    otu_model = micom.util.join_models(value, str(key))
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/Hale_2018_model_building_and_combining/Hale_2018_otu_specific_models/")
    cobra.io.write_sbml_model(otu_model, filename)
    del otu_model
    counter= counter+1
    print(str(counter))
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/AGORA_models/")
del counter, filename, key, value

#ugh. not working on MSI. SPlit dataframe

agora_df_1 = all_agora_merged_condensed2.iloc[:172,:]
agora_df_2 = all_agora_merged_condensed2.iloc[172:344,:]
agora_df_3 = all_agora_merged_condensed2.iloc[344:516,:]
agora_df_4 = all_agora_merged_condensed2.iloc[516:688,:]
agora_df_5 = all_agora_merged_condensed2.iloc[688:860,:]
agora_df_6 = all_agora_merged_condensed2.iloc[860:1032,:]
agora_df_7 = all_agora_merged_condensed2.iloc[1032:1204,:]
agora_df_8 = all_agora_merged_condensed2.iloc[1204:,:]

os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Hale_2018_model_picking/")
agora_df_1.to_pickle('AGORA_models_that_need_OTU_builds_set1.pickle')
agora_df_2.to_pickle('AGORA_models_that_need_OTU_builds_set2.pickle')
agora_df_3.to_pickle('AGORA_models_that_need_OTU_builds_set3.pickle')
agora_df_4.to_pickle('AGORA_models_that_need_OTU_builds_set4.pickle')
agora_df_5.to_pickle('AGORA_models_that_need_OTU_builds_set5.pickle')
agora_df_6.to_pickle('AGORA_models_that_need_OTU_builds_set6.pickle')
agora_df_7.to_pickle('AGORA_models_that_need_OTU_builds_set7.pickle')
agora_df_8.to_pickle('AGORA_models_that_need_OTU_builds_set8.pickle')

#do 5-8 if this works
#do the same thing with mambo, to create some consistency in model gen
mambo_otus = pd.read_pickle("MAMBO_models_that_need_OTU_builds.pickle")
mambo_df_1 = mambo_otus.iloc[:32,:]
mambo_df_2 = mambo_otus.iloc[32:,:]
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Hale_2018_model_picking/")
mambo_df_1.to_pickle('MAMBO_models_that_need_OTU_builds_set1.pickle')
mambo_df_2.to_pickle('MAMBO_models_that_need_OTU_builds_set2.pickle')


carveme_otus = pd.read_pickle("CarveMe_models_that_need_OTU_builds.pickle")

########keeping track of all OTU models generated.
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/Hale_2018_model_building_and_combining/Hale_2018_otu_specific_models/")

for file in os.listdir():
    print(file)
    

os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Hale_2018_model_picking/")
carveme_otus = pd.read_pickle("CarveMe_models_that_need_OTU_builds.pickle")
agora_otus = pd.read_pickle("AGORA_models_that_need_OTU_builds.pickle")

##########this is nothign, just testing something
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Building_MICOM_communities_Hale2018/")
my_file = open("Hale2018_OTU_abundances_2.txt", "r")
content_list = my_file.readlines()
print(content_list)


















