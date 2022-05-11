#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 16:39:07 2020

@author: adamo010
"""

#The purpose of this code is to generate .pickle MICOM community files for each patient in
#the Burns dataset (88 communities total, 2 sites per 44 patients).

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

#I realize that I'll need to re-do abundance filtering based on the OTU table that I actually ended up with
#bring in absolute abundance OTU table and do relative abundances from a filtered version of that.

#importing files
model_ids = pd.read_csv("06.22.20_summary_of_otus_with_models.csv")
preprocessed_otus = pd.read_csv("Contaminant_removal/Burns_absolute_OTU_table.csv") #note that this is in /Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal

preprocessed_otus.rename(columns={"#OTU_ID":"OTU_ID"}, inplace = True)

#create a list of OTU_IDs in model_IDs
OTUs_with_models = []
for row in model_ids["OTU_ID"]:
    OTUs_with_models.append(row)
del row

#create a new dataframe where only OTU_IDs matching the list of OTU_IDs in model_IDs is included
filtered_absolute_abundance_table = preprocessed_otus[preprocessed_otus['OTU_ID'].isin(OTUs_with_models)]

#convert the table to relative abundances.
#step 0: pull all columns that need to be summed (i.e. all OTU columns)- will use a lot
sample_names = []            #this is a list of all the column names that will need to be scanned for values >0.001
for col in filtered_absolute_abundance_table.columns:       #shouldn't matter which dataframe we use here but use raw one just in case
    if 'Sample_' in str(col):
        sample_names.append(str(col))
    del col

#step 1: collapse by class- not really applicable here but the code won't run without it
otu_table_collapsed = filtered_absolute_abundance_table.groupby(['taxonomy', 'OTU_ID'])[sample_names].sum()
#fine, but we actually need everything else in taxsplit_otus that's not in sample_names
#or IS THERE???
#maybe we just keep taxonomy_orig and #OTU_ID. the rest of the taxonomy stuff we can add again later.

#step 2: calculate relative abundances
otu_table_relabund = otu_table_collapsed .loc[:,].div(otu_table_collapsed .sum(axis=0))
#it is worth noting that this works BECAUSE in the groupby function in step 1, taxonomy_orig and #OTU_ID are both
#set as the index. Otherwise, can't do maths on strings 

#however, you DO have to reset the index for downstream stuff.
otu_table_relabund.reset_index(inplace=True)

#ok, so now OTU_table_relabund has all relative abundances.
#now, I need to attach model file names to each otu
#probably the best way to do that is create a dictionary of keys(otus) and values (model file names)
#have I done this before? Yes. Do I remember how? No.

model_file_name_dict = pd.Series(model_ids.final_model_name.values,index=model_ids.OTU_ID).to_dict()
#nice

#use dictionary to create a new model_file_name column
otu_table_relabund['model_file_name'] = otu_table_relabund['OTU_ID'].map(model_file_name_dict)

#save this guy as a csv JUST IN CASE we need to look at it again. Which we definitely will.
otu_table_relabund.to_csv("Individual_specific_microbiomes/08.13.20_relabund_OTU_table_for_OTUs_with_models.csv")

#all right, now for the real hard stuff. Create several new data frames, one for each column, where only non-zero abundance rows are included.
#above_35 = titanic[titanic["Age"] > 35]
#sample_name_df = otu_table_relabund[otu_table_relabund["sample_XX"] > 0]

#the things I want to do are:
#1) create new dataframes, one for each sample
#2) for each dataframe, only keep non-zero abundance rows
#3) for each dataframe, only keep the following columns: taxonomy, OTU_id, model_file_name, Sample_XX(relative abundance))

#probably have to create a function for this.
#I KNOW you're not supposed to iteratively create variables, but I'm doing it anyway. 

def subsetting_OTU_table(sample_name):
    subset_df = otu_table_relabund[otu_table_relabund[str(sample_name)] > 0]
    subset_df2 = subset_df[["taxonomy", "OTU_ID", "model_file_name", str(sample_name)]]
    subset_df2.to_csv("Individual_specific_microbiomes/_abundances_and_models.csv")
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes")
    for file in os.listdir():                   
        src=file
        if fnmatch.fnmatch(file, "_abundances_and_models.csv"):
            dst = str(sample_name)+file
            os.rename(src,dst)
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/")       
    return

for elem in sample_names:
    subsetting_OTU_table(str(elem))

#so this works, but something is fishy. All output files are either 874 bytes or 4kb. The smaller files have some taxa
    #listed but also some are just dots.

#yeah, ok, something is really wrong here. Something in the transition of the subsetted OTU table to a csv    
#work on it here with Sample_27 
#subset_df = otu_table_relabund[otu_table_relabund[str("Sample_27")] > 0]
#subset_df2 = subset_df[["taxonomy", "OTU_ID", "model_file_name", str("Sample_27")]]
#subset_df2.to_csv("Individual_specific_microbiomes/Dummy.csv")
#output = open("Individual_specific_microbiomes/Dummy.csv", "w")
#output.write(str(subset_df2))
#output.close()

#So I'm an ass and I'm not going to delete this because I think it's an important learning moment.
#this whole code below
    #output = open("Individual_specific_microbiomes/Dummy.csv", "w")
    #output.write(str(subset_df2))
    #output.close()
#should be written as this:
    #subset_df2.to_csv("Individual_specific_microbiomes/Dummy.csv")
#I don't know what the hell I was thinking in that old microbiome_fba_on_diets.py file, but shit.

#OOOOOOKAY, moving on.
#we now have a set of tables containing the necessary information for MICOM. These files are all named as follows:
#Sample_XX_abundances_and_models.csv
#Now we will import these and make communities with them (and very probably break my poor computer).

#let's start with Sample 27.
Sample_27_df = pd.read_csv("Individual_specific_microbiomes/Sample_27_abundances_and_models.csv")    
Sample_27_df.drop(["Unnamed: 0", "taxonomy"], axis=1, inplace=True)    #get rid of extra column.
    
#see MICOM_healthy_vs_CRC_trial.py for info on these next steps
    
from micom import Community    
#sample_27_community = Community(Sample_27_df, solver='gurobi')    
#can't do this raw; need an id column
Sample_27_df.rename(columns={'OTU_ID': 'id', 'Sample_27': 'abundance', 'model_file_name': 'file'}, inplace=True)
#I keep getting this "AttributeError: Can only use .str accessor with string values!" error. Maybe the stupid
#program is reading my id column as not strings?
#df['DataFrame Column'] = df['DataFrame Column'].astype(float)
Sample_27_df['id']= Sample_27_df["id"].astype(str)
#so, that worked, but now need to change directory to /Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/OTU_models
#there is probably a way to do this in the command, but I don't know it. 
sample_27_community = Community(Sample_27_df, solver='gurobi')     
#and now we wait. Looking at 10-20 seconds per model, this particular sample has 58 OTUs, so... probably about 20 minutes.
#yep. About 22 minutes.
#now, save to pickle.
sample_27_community.to_pickle("sample_27_community.pickle")
#cool.

#in the meantime, I'll work on writing parallelization code.
#first, create a list of all the files with each sample's OTU table    
sample_OTU_table_files = []
for file in glob.glob("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/*abundances_and_models.csv"):
    sample_OTU_table_files.append(file)   
del file
sample_OTU_table_samplenames = []
for elem in sample_OTU_table_files:
    elem = elem.lstrip("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/")
    elem = elem.rstrip("_abundances_and_models.csv")
    sample_OTU_table_samplenames.append(elem)
del elem

#merge these into a dictonary; then we can iterate over both sets.
sample_dict = {sample_OTU_table_samplenames[i]: sample_OTU_table_files[i] for i in range(len(sample_OTU_table_samplenames))}

def creating_MICOM_tables(sample_file):
    file_name = str(sample_file)
    sample_name = file_name.lstrip("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/").rstrip("_abundances_and_models.csv")
    sample_name.df = pd.read_csv("Individual_specific_microbiomes/sample_file")

#step 1: importing and editing all dataframes
edited_sample_dict = {}
for key, value in sample_dict.items():
    dataframe = pd.read_csv(value)
    dataframe.rename(columns = {'OTU_ID': 'id', key: 'abundance', 'model_file_name': 'file'}, inplace=True)
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
#seems like we have 76 pickled communities and one unnamed pickled community (probably where it got stuck).
#Yep, sure enough. See 08.21.20_analysis_of_MICOM_community_building_time.xls Pickling tab. Based on the order that the samples were
#listed, Sample 61 is where things timed out. (note that they are not in numerical order)        

#I'm going to try to make a list of the missing samples, subset the dictionary finished_communities by that list, and rerun the to_pickle code.
        
missing_samples = ["Sample_61","Sample_66", "Sample_11", "Sample_48", "Sample_88", "Sample_19", "Sample_80","Sample_37","Sample_40","Sample_95","Sample_22","Sample_73"]        

#expectedResult = [d for d in exampleSet if d['type'] in keyValList]
#res = [dictA[i] for i in key_list if i in dictA]
round2_finished_communities = {key: value for key, value in finished_communities.items() if key in missing_samples}

#okay, so we're at the point where none of the variables are working. Time to close down Spyder, reboot, and start this ALL OVER AGAIN.
#did everything up to step 1 above.

#now, filter edited_sample_dict by missing_samples so we're not rerunning MICOM on everything.
round2_edited_sample_dict = {key: value for key, value in edited_sample_dict.items() if key in missing_samples}

#NOW, redo step 2. Started 
finished_communities_round2 = {}
for key, value in round2_edited_sample_dict.items():
    sample_community = Community(value, solver='gurobi')
    finished_communities_round2.update( {key: sample_community})     

#redoing step 3    
for key, value in finished_communities_round2.items():
    value.to_pickle("_community.pickle")
    for file in os.listdir():
        src= file
        if fnmatch.fnmatch(file, "_community.pickle"):
            dst = str(key)+file
            os.rename(src, dst)    

    