#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 11:04:43 2020

@author: adamo010
"""
#updates from 05.28.20: have now tracked down all otus with models


#So, now that we have merged/Carved metabolic models for a bunch of taxa (plus a couple more contaminants), time to go through
#and create a 'slimmed' relative abundance table that only has OTUs with models in it.

#this code is somewhat poached from 04.27.20_OTU_table_filtering_pipeline_clean.py


import scipy
import numpy as np
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
import ast

#I think the best course of action is the following:
#1) make a list of all the otus with models. Scrape this information from the names of the .xml files with otu names
#2) update the contaminants list to include the OTUs ided as contaminants in 05.01.20_OTUs_to_build_models_for.csv
#3) actually create an OTU table for contaminants (rather than a list of otus or indices)
#4) get all the OTUs to remove in order with reasons on why they're being removed
#5) get metadata on OTUs with models
#6) filter the absolute abundance table to something that only includes the otus with models- see 05.20.20_contaminant_otu_table.py
#7) calculate relative abundances for model-having otus
#8) compare the model-having to model-absent otus.
#then, further down the line:
#9) identify samples that are CRC vs non-CRC
#10) figure out if there is any bias in model-having vs model-absent OTUs in CRC vs non-CRC samples. 

#Step 1. make a list of all the model files for specific OTUs
from os import listdir
from os.path import isfile, join
onlyfiles1 = [f for f in listdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining") if isfile(join("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining", f))]
onlyfiles2 = [f for f in listdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/0_genome_seqs_for_model_building") if isfile(join("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/0_genome_seqs_for_model_building", f))]

otu_model_list = []
for elem in onlyfiles1:
    if elem.endswith("model.xml"):
        otu_model_list.append(elem)
del elem    

for elem in onlyfiles2:
    if elem.endswith("model.xml"):
        otu_model_list.append(elem)
del elem     

#extract a list of the otus which have models
otus_with_models = []
for elem in otu_model_list:
    elem = elem.rstrip("model.xml")
    elem = elem.lstrip("otu_")
    otus_with_models.append(elem)
del elem    

#step 2-3: update the contaminant models list and make an otu table
#oh boy, this is going to take some time. do this in 05.28.20_contaminant_otu_table.py
#at last, on 06.03.20, I have that information. Yikes alive, that took a while. In my defense, it was really hard code and also
#George Floyd got murdered by cops on 05.25.20
contaminant_otus = pd.read_csv("Contaminant_removal/06.16.20_table_of_contaminant_OTUs.csv")

#import some files I'll need to do my caluclations.
preprocessed_otus = pd.read_csv("Contaminant_removal/Burns_absolute_OTU_table.csv")
otus_with_no_genomes = pd.read_csv("Model_building_and_combining/05.28.20_otus_with_no_genome_seqs.csv")
#contaminant_otus = pd.read_csv("Contaminant_removal/05.28.20_Burns_contaminant_list.txt", sep="    ", header = None, engine='python')
otus_with_carved_models = pd.read_csv("Model_building_and_combining/05.01.20_OTUs_to_build_models_for.csv")
otus_with_agora_models1 = pd.read_csv("Model_picking/06.16.20_OTUs_with_models_V1.csv")
otus_with_agora_models2 = pd.read_csv("Model_picking/05.01.20_OTUs_with_models_V1_taxonomy_fixed.csv")
otus_with_mambo_models = pd.read_csv("Model_picking/06.16.20_OTUs_with_mambo_models_V1_plus2.csv")   

#otus_without_models will consist of otus_with_no_genomes, contaminant_otus

#step 4a: metadata for otus without models
#do this in dictionary form: keys have to be unique, so they should be the otu_ids
otus_with_no_genomes.rename(columns = {"#OTU_ID": "OTU_ID"}, inplace=True)
otu_keys = []
reason_left_out_values = []
for row in otus_with_no_genomes["OTU_ID"]:
    otu_keys.append(row)
    reason_left_out_values.append("no genome or model available")
del row

no_genomes_dict = dict(zip(otu_keys, reason_left_out_values)) 
del otu_keys
del reason_left_out_values

#step 4b: metadata for contaminant otus
#done in step 2-3
#now the trick is to do the same thing as no_genomes_dict: just a series of key-value pairs
otu_keys = []
reason_left_out_values = []
for row in contaminant_otus["OTU_ID"]:
    otu_keys.append(row)
    reason_left_out_values.append("contaminant")
del row

contaminants_dict = dict(zip(otu_keys, reason_left_out_values)) 
del otu_keys
del reason_left_out_values

#step 5. Get information on the taxonomic level and number of models that went into each otu model.
#this will probably be different for every input file
#goal: for each of the last four series, get a pandas series with the following columns: otu_id, model_source, num_source_models, final_model_name 
#goal2: for each of the first two series, get a pandas series with the following columns: otu_id, reason_left_out

#otus_with_models will consit of otus_with_carved_models, otus_with_agora_models1, otus_with_agora_models2, otus_with_mambo_models
otus_with_carved_models.rename(columns = {"#OTU_ID": "OTU_ID"}, inplace=True)
otus_with_carved_models.rename(columns = {"File name": "filename"}, inplace=True)
otus_with_agora_models1.rename(columns = {"#OTU_ID": "OTU_ID"}, inplace=True)
otus_with_agora_models2.rename(columns = {"#OTU_ID": "OTU_ID"}, inplace=True)
otus_with_mambo_models.rename(columns = {"#OTU_ID": "OTU_ID"}, inplace=True)

print(otus_with_carved_models.columns)
#for this one, need to count the number of elements in the column File name
#drop rows where no genome files were found (so no models were made): should eliminate 5 columns
otus_with_carved_models.dropna(subset = ['filename'], inplace = True)

otu_keys = []
number_of_genomes = []
for row in otus_with_carved_models["OTU_ID"]:
        otu_keys.append(row)

#okay, need to turn File name column into a list
otus_with_carved_models['filename'] = otus_with_carved_models.filename.apply(lambda x: x.strip().split(';'))
        
for row in otus_with_carved_models["filename"]:
    count = 0
    for elem in row:
        count = count+1
    #print(count)
    number_of_genomes.append(count)
del row
del count

carved_models_dict = dict(zip(otu_keys, number_of_genomes)) 
del number_of_genomes
#nice.

#now, I can get the final model name off of the otu_id column

#let's start our dataframe-to-be-merged
#use otu_keys as constant keys for a bunch of dictionaries
model_source_values = []
final_model_name_values = []
for elem in otu_keys:
    model_source_values.append("CarveMe")
    final_model_name_values.append("otu_"+str(elem)+"model.xml")
del elem    

carved_model_source_dict = dict(zip(otu_keys, model_source_values)) 
carved_model_name_dict = dict(zip(otu_keys, final_model_name_values))
del model_source_values
del final_model_name_values

#create our data frame!
carved_models_summary = pd.DataFrame(list(carved_models_dict.items()),columns = ['OTU_ID','number_of_models_combined']) 
carved_models_summary["model_source"]= carved_models_summary['OTU_ID'].map(carved_model_source_dict)
carved_models_summary["final_model_name"]= carved_models_summary['OTU_ID'].map(carved_model_name_dict)

#oh crap, I forgot about getting the taxonomic level.
#I'm just going to have to edit the input file and go from there
print(otus_with_carved_models.columns)
otus_with_carved_models.rename(columns = {"Level at which to make the model": "model_phylogeny_level"}, inplace=True)

carved_model_phylogeny_values = []
for row in otus_with_carved_models["model_phylogeny_level"]:
    carved_model_phylogeny_values.append(str(row))
del row

carved_model_phylogeny_dict = dict(zip(otu_keys, carved_model_phylogeny_values)) 
carved_models_summary["model_taxonomic_level"]= carved_models_summary['OTU_ID'].map(carved_model_phylogeny_dict)

##########COOL! that works.
print(carved_models_summary.columns)
#need these columns in this order moving forward: 'OTU_ID', 'number_of_models_combined', 'model_source','final_model_name', "model_taxonomic_level')

#next up: mambo models
#borrow some code from 05.01.20_combining_models_for_Burns_OTUs
print(otus_with_mambo_models.columns)
otus_with_mambo_models.drop(["Unnamed: 0", "Unnamed: 0.1", "Unnamed: 0.1.1"], axis = 1, inplace=True) #drop 4 boring columns
#create a new 'models list' column and merge all data from other model columns
otus_with_mambo_models["model_columns_list"] = otus_with_mambo_models["species_matched_model_list"].astype(str) +otus_with_mambo_models["genus_matched_model_list"].astype(str) +otus_with_mambo_models["family_matched_model_list"].astype(str)
#that... kind of worked? Get rid of nan's that were appended to the list
otus_with_mambo_models["model_columns_list"] = otus_with_mambo_models["model_columns_list"].map(lambda x: x.lstrip('nan').rstrip('nan'))

#now, create a column where the phylogenic level of the model is identified
#first, how do we figure out the level at which a model is identified?
#I have three columns (species/genus/family_matched_model_list). Depending on which one is NOT blank, I want to return
#a specific value and append that value to a new column corresponding to the appropriate row.
#the problem is that these rows contain (float) nans and (string) model names; no function I try can handle both.
#what I would like is 'if species is string, print "species"

#fuck it. Convert everything to strings first
otus_with_mambo_models['species_matched_model_list']= otus_with_mambo_models['species_matched_model_list'].astype(str)
otus_with_mambo_models['genus_matched_model_list']= otus_with_mambo_models['genus_matched_model_list'].astype(str)
otus_with_mambo_models['family_matched_model_list']= otus_with_mambo_models['family_matched_model_list'].astype(str)

def phylogeny_level_test(row):
  if row['species_matched_model_list'] != "nan":
    return "Species"
  elif row['genus_matched_model_list'] != "nan":
    return "Genus"
  elif row['family_matched_model_list'] != "nan":
    return "Family"

otus_with_mambo_models['model_taxonomic_level'] = otus_with_mambo_models.apply(lambda row: phylogeny_level_test(row), axis=1)

#otus_with_mambo_models.drop(["model_taxonomic_level"], axis = 1, inplace=True) 
    
#genuinely, this took me three fucking days to figure out. 

#######okay, now we create a column where we have counts of the models
#need to turn File name column into a list
otus_with_mambo_models['model_columns_list'] = otus_with_mambo_models['model_columns_list'].apply(ast.literal_eval)

for row in otus_with_mambo_models["model_columns_list"]:
    print(type(row))

#otus_with_mambo_models2 = copy.deepcopy(otus_with_mambo_models)

def counting_models(row):
    counter = 0
    for elem in row["model_columns_list"]:
        counter= counter+1
    return counter    

otus_with_mambo_models['number_of_models_combined'] = otus_with_mambo_models.apply(lambda row: counting_models(row), axis=1)

#cool.

#add a model source column
def model_source(row):
    return "MAMBO"

#otus_with_mambo_models2 = copy.deepcopy(otus_with_mambo_models)
otus_with_mambo_models['model_source'] = otus_with_mambo_models.apply(lambda row: model_source(row), axis=1)

#this whole function-apply thing is slick once it works.
#all right, just need final_model_name now.
#probably do this as a list because... it's easier.

model_name_list = []
for elem in otus_with_mambo_models["OTU_ID"]:
    model_name_list.append("otu_"+str(elem)+"model.xml")
del elem
  
#otus_with_mambo_models2 = copy.deepcopy(otus_with_mambo_models)
otus_with_mambo_models['final_model_name'] = model_name_list
del model_name_list

#all right, now we have OTU_ID, model_taxonomic_level, number_of_models_combined, model_source, and final_model_name
#need to export these to a new dataframe called mambo_models summary.
#need these columns in this order moving forward: 'OTU_ID', 'number_of_models_combined', 'model_source','final_model_name', "model_taxonomic_level')

mambo_models_summary = otus_with_mambo_models[['OTU_ID', 'number_of_models_combined', 'model_source','final_model_name', 'model_taxonomic_level']].copy()

#killin it.
#okay, now to agora 1 and agora2


############Agora1#############
#step 1: make a total models column
print(otus_with_agora_models1.columns)
otus_with_agora_models1.drop(["Unnamed: 0", "Unnamed: 0.1"], axis = 1, inplace=True) #drop 4 boring columns
#create a new 'models list' column and merge all data from other model columns
otus_with_agora_models1["model_columns_list"] = otus_with_agora_models1["species_matched_model_list"].astype(str) +otus_with_agora_models1["genus_matched_model_list"].astype(str) +otus_with_agora_models1["family_matched_model_list"].astype(str)
#that... kind of worked? Get rid of nan's that were appended to the list
otus_with_agora_models1["model_columns_list"] = otus_with_agora_models1["model_columns_list"].map(lambda x: x.lstrip('nan').rstrip('nan'))

#step 2 from HELL: create a column where the phylogenic level of the model is identified

#fuck it. Convert everything to strings first
otus_with_agora_models1['species_matched_model_list']= otus_with_agora_models1['species_matched_model_list'].astype(str)
otus_with_agora_models1['genus_matched_model_list']= otus_with_agora_models1['genus_matched_model_list'].astype(str)
otus_with_agora_models1['family_matched_model_list']= otus_with_agora_models1['family_matched_model_list'].astype(str)

#function phylogeny_level_test is defined in mambo stuff
otus_with_agora_models1['model_taxonomic_level'] = otus_with_agora_models1.apply(lambda row: phylogeny_level_test(row), axis=1)

#step 3: create a column where we have counts of the models
otus_with_agora_models1['model_columns_list'] = otus_with_agora_models1['model_columns_list'].apply(ast.literal_eval)

#function counting_models defined above in mambo stuff
otus_with_agora_models1['number_of_models_combined'] = otus_with_agora_models1.apply(lambda row: counting_models(row), axis=1)

#step 4: add a model source column
#function model_source HAS TO BE RE-DEFINED HERE
def model_source(row):
    return "AGORA"

otus_with_agora_models1['model_source'] = otus_with_agora_models1.apply(lambda row: model_source(row), axis=1)

#step 5: create a model name list and append it as a column

model_name_list = []
for elem in otus_with_agora_models1["OTU_ID"]:
    model_name_list.append("otu_"+str(elem)+"model.xml")
del elem
  
otus_with_agora_models1['final_model_name'] = model_name_list
del model_name_list

#step 6: make a summary dataframe
agora1_models_summary = otus_with_agora_models1[['OTU_ID', 'number_of_models_combined', 'model_source','final_model_name', 'model_taxonomic_level']].copy()

#goddamn it. Nailed it except for the model_source function, which I stupidly left as MAMBO  . At least it was easy to fix.   
    
############Agora2#############
#step 1: make a total models column
print(otus_with_agora_models2.columns)
otus_with_agora_models2.drop(["Unnamed: 0", "Unnamed: 0.1"], axis = 1, inplace=True) #drop 4 boring columns
#create a new 'models list' column and merge all data from other model columns
otus_with_agora_models2["model_columns_list"] = otus_with_agora_models2["species_matched_model_list"].astype(str) +otus_with_agora_models2["genus_matched_model_list"].astype(str) +otus_with_agora_models2["family_matched_model_list"].astype(str)
#that... kind of worked? Get rid of nan's that were appended to the list
otus_with_agora_models2["model_columns_list"] = otus_with_agora_models2["model_columns_list"].map(lambda x: x.lstrip('nan').rstrip('nan'))

#step 2 from HELL: create a column where the phylogenic level of the model is identified

#fuck it. Convert everything to strings first
otus_with_agora_models2['species_matched_model_list']= otus_with_agora_models2['species_matched_model_list'].astype(str)
otus_with_agora_models2['genus_matched_model_list']= otus_with_agora_models2['genus_matched_model_list'].astype(str)
otus_with_agora_models2['family_matched_model_list']= otus_with_agora_models2['family_matched_model_list'].astype(str)

#function phylogeny_level_test is defined in mambo stuff
otus_with_agora_models2['model_taxonomic_level'] = otus_with_agora_models2.apply(lambda row: phylogeny_level_test(row), axis=1)

#step 3: create a column where we have counts of the models
otus_with_agora_models2['model_columns_list'] = otus_with_agora_models2['model_columns_list'].apply(ast.literal_eval)

#function counting_models defined above in mambo stuff
otus_with_agora_models2['number_of_models_combined'] = otus_with_agora_models2.apply(lambda row: counting_models(row), axis=1)

#step 4: add a model source column
#function model_source HAS TO BE RE-DEFINED HERE
def model_source(row):
    return "AGORA"

otus_with_agora_models2['model_source'] = otus_with_agora_models2.apply(lambda row: model_source(row), axis=1)

#step 5: create a model name list and append it as a column

model_name_list = []
for elem in otus_with_agora_models2["OTU_ID"]:
    model_name_list.append("otu_"+str(elem)+"model.xml")
del elem
  
otus_with_agora_models2['final_model_name'] = model_name_list
del model_name_list

#step 6: make a summary dataframe
agora2_models_summary = otus_with_agora_models2[['OTU_ID', 'number_of_models_combined', 'model_source','final_model_name', 'model_taxonomic_level']].copy()

#it is genuinely a real mixed bag when the code you've been working on for a week suddenly becomes something that takes 10 seconds to replicate

#####################OOOOOOOOOOOOOOKKKKKKKKKKAAAAAAAAAAAAYYYYYYYYYYYYYY time to smash some summary dataframes together.

#otus with models:
frames_with_models = [agora1_models_summary, agora2_models_summary, carved_models_summary, mambo_models_summary]
summary_of_all_otus_with_models= pd.concat(frames_with_models)

summary_of_all_otus_with_models.to_csv("06.23.20_summary_of_otus_with_models.csv")
   
#otus_without_models:
no_genomes_dict  
contaminants_dict  

from pandas import DataFrame
no_genomes_df = DataFrame(list(no_genomes_dict.items()),columns = ['OTU_ID','reason_for_no_model']) 
contaminants_df = DataFrame(list(contaminants_dict.items()),columns = ['OTU_ID','reason_for_no_model']) 

frames_with_models2 = [no_genomes_df, contaminants_df]
summary_of_all_otus_without_models= pd.concat(frames_with_models2)

summary_of_all_otus_without_models.to_csv("06.23.20_summary_of_otus_without_models.csv")









