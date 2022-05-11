#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 16:42:41 2020

@author: adamo010
"""
#this is based on 05.01.20_combining_models_for_Burns_otus.py; I am currently using models that I made from CarveMe.

#edits from 04.29: added a few more OTUs to the model-matched database, so had to include them here.
#also, adding ALL OTUs with matched models (not just the ones with multiple models per OTU), since it turns out that
#there is (surprise, surprise) a tonne of database management to do before pulling models. 


#the goal of this code is to a) make a list of all the AGORA models to download, since they don't all download at once,
#b) combine all the models together where multiple models have been assigned to a single OTU
#c) put all models for this particular dataset in a central repository for running later.

#problems I forsee: mostly filepath ones. Honestly. You know, stupid problems.

import micom
from micom.util import (
    load_model,
    join_models,
    add_var_from_expression,
    adjust_solver_config,
    clean_ids,
)
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
import ast

#the tricky bit will be tying each combined model to the correct OTU.
#I think the best approach is to collapse the species/genus/family_level_models columns into a single column. Not sure how to deal with 
#AGORA vs MAMBO models, though. Shouldn't be an issue?

otus_with_models = pd.read_csv("05.01.20_OTUs_to_build_models_for.csv")
#clean up this file:
#otus_with_models.drop(["Unnamed: 0", "Unnamed: 0.1", "Unnamed: 0.1.1", "Unnamed: 0.1.1.1"], axis = 1, inplace=True) #drop 4 boring columns
otus_with_models.rename(columns = {"#OTU_ID": "OTU_ID"}, inplace=True)
otus_with_models.rename(columns = {"File name": "model_files"}, inplace=True)

#delete rows with nans in list_of_model_files column 
otus_with_models.dropna(subset = ["model_files"], inplace=True)

#split the model_files column contents into lists
model_list = []
for row in otus_with_models["model_files"]:
    row= str(row)
    print(type(row))
    row= row.split(sep=';')
    model_list.append(row)
del row    

otus_with_models['list_of_models'] = model_list

#create a new model_name column based on model_files where the _protein.faa is replaced with .xml
model_file_list = []

for row in otus_with_models["list_of_models"]:
    row_model_list = []
    while(" " in row):  #this removes blank line that screwed up some stuff later on.
        row.remove(" ")
    for elem in row:
        #print(elem)
        elem.strip()
        elem = re.sub('_protein.faa$', '.xml', elem)
        elem = re.sub('_genomic.fna$', '.xml', elem)
        elem = elem.strip()
        row_model_list.append(elem)
    model_file_list.append(row_model_list)    

otus_with_models['list_of_model_files']= model_file_list

#make a list of all the model files
from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/0_genome_seqs_for_model_building") if isfile(join("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/0_genome_seqs_for_model_building", f))]

modelfiles = []
for elem in onlyfiles:
    if elem.endswith(".xml"):
        modelfiles.append(elem)

modelfiles_to_combine = []

for row in otus_with_models['list_of_model_files']:
    for elem in row:
        print(elem)
        if elem not in modelfiles_to_combine:
            modelfiles_to_combine.append(elem)

#okay, so it turns out we do have to do that fun name-matching thing. 
missing_models = np.setdiff1d(modelfiles_to_combine, modelfiles)
# yields the elements in `list_2` that are NOT in `list_1`

#only two are fucked up, which is good.
#oh right, the Paraoerskovia thing had a missing underscore that I changed in the filename but not in the table. 
#Rubritalea profundi doesn't have any underscores in the modelfiles. 
#Normally I would just change the filenames, but I've already been editing the input file, so I'm just going to do that again. 
#if this is working, missing_models should be empty

#okay, one of the rows in my list_of_model_files list is blank. need to fix that. 


#the below code should work but running into micom issues. 
#this will theoretically run join_models.    
#remember to switch to the correct directory here

otu_list = []               #create an OTU list to serve as keys for model file list dictionary
for row in otus_with_models['OTU_ID']:
    otu_list.append(row)
del row

combined_model_list = []

for row in otus_with_models['list_of_model_files']:
    print(row)
    #print("      ")
    #print(model_list)
    otu_model = micom.util.join_models(row, "otu_model")
    combined_model_list.append(otu_model)

model_files_dict = dict(zip(otu_list, combined_model_list)) #good thing lists are ordered, that's all I'll say. 

#so, now the trick is to iterate through the dictionary and save each model as an xml file. 
for key, value in model_files_dict.items():
    filename = "otu_"+str(key)+"model.xml"
    cobra.io.write_sbml_model(value, filename)  #I don't know if this will work, but maybe?

########Borrow from 05.22.20_testing_CarveMe_and_combo_models to check these
#get all merged models into a list called merged_models
onlyfiles2 = [f for f in listdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/0_genome_seqs_for_model_building") if isfile(join("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/0_genome_seqs_for_model_building", f))]

merged_models = []
for filename in onlyfiles2:
    if filename.endswith("model.xml"):
        merged_models.append(filename)
del filename
del onlyfiles2

#cool. Now, I want to run FBA on all the models from carved_models and merged_models
#borrow some of the model-running code from microbiome_fba_cleancopy.py

model_test_results = {}

def metabolic_model(taxon_model):
    fba_model= cobra.io.read_sbml_model(taxon_model)
    fba_model_sol = fba_model.optimize()
    sol_string = str(fba_model_sol)
    model_test_results.update( {taxon_model : sol_string} )
    return

for model in merged_models:
    metabolic_model(model)
    
#save models and take a look at the results
merged_model_results = open("05.27.20_carved_and_merged_model_test_results.csv", "w")
writer= csv.writer(merged_model_results)
for key, value in model_test_results.items():
    writer.writerow([key, value])
    
merged_model_results.close()    

#Nice! All those models ran too.
