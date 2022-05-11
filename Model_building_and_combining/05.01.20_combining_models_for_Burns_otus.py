#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 16:42:41 2020

@author: adamo010
"""
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

otus_with_models = pd.read_csv("05.01.20_OTUs_that_need_merged_models.csv")
#clean up this file:
otus_with_models.drop(["Unnamed: 0", "Unnamed: 0.1", "Unnamed: 0.1.1", "Unnamed: 0.1.1.1"], axis = 1, inplace=True) #drop 4 boring columns
otus_with_models.rename(columns = {"#OTU_ID": "OTU_ID"}, inplace=True)

#create a new 'models list' column and merge all data from other model columns
otus_with_models["model_columns_list"] = otus_with_models["species_matched_model_list"].astype(str) +otus_with_models["genus_matched_model_list"].astype(str) +otus_with_models["family_matched_model_list"].astype(str)

#that... kind of worked? Get rid of nan's that were appended to the list
otus_with_models["model_columns_list"] = otus_with_models["model_columns_list"].map(lambda x: x.lstrip('nan').rstrip('nan'))
#AND get rid of existing.xmls?

#next steps: split model_columns_list into indiv models (vs lists of models)
otus_with_models['model_columns_list'] = otus_with_models['model_columns_list'].apply(ast.literal_eval)
#the ast module is new for me, but it looks like it worked
for row in otus_with_models['model_columns_list']:
    print(type(row))
del row    
#now we have lists instead of strings
    
#FOR FUCK'S SAKE the model names I downloaded all have underscores in them. also need to remove the equals signs
cleaned_otus_with_models = copy.deepcopy(otus_with_models)    

#update: some of the models (basically the MAMBO models) now have .xml.xml at the end of their filenames. Which is great
#I tried to if/then into in the above iteration, but it didn't work and I don't want to fuck with this too much.
#try removing .xml extensions from filenames instead of skipping them in the below 

otu_list = []               #create an OTU list to serve as keys for model file list dictionary
for row in cleaned_otus_with_models['OTU_ID']:
    otu_list.append(row)
del row

list_of_file_lists = []     #create an list of model file lists to serve as values for model file list dictionary
for row in cleaned_otus_with_models['model_columns_list']:
    file_list = []      #initiate an empty list: for each row, the model files will be added to this list
    for elem in row:        #first, have to clean up the model names so they match file names
        elem_temp = elem.split(" =",1)
        elem_temp2 = elem_temp[0]
        elem_temp3 = elem_temp2.strip() #take off white space at the ends
        elem_temp4 = elem_temp3.rstrip(".xml")      #have to pull .xml off MAMBO models
        elem_temp5 = elem_temp4.replace(", ","_").replace(". ","_").replace(".","_").replace(" ","_").replace("-","_").replace("/","_").replace(":","_").replace("(","").replace(")","").replace("__","_").replace("'","")
        elem_temp6 = elem_temp5+ ".xml"
        #print(elem_temp4)
        file_list.append(elem_temp6)
        #print(file_list)    
    list_of_file_lists.append(file_list)    #append the list of files to the list of lists of files (what has my life come to?)
del elem
del elem_temp
del elem_temp2
del elem_temp3
del elem_temp4
del elem_temp5
del elem_temp6
del file_list
del row    


#make a dictionary where otus are matched with their lists of associated model files
model_files_dict = dict(zip(otu_list, list_of_file_lists)) #good thing lists are ordered, that's all I'll say. 

#dunk this as a new column in cleaned_otus_with_models
cleaned_otus_with_models['model_file_list'] = cleaned_otus_with_models['OTU_ID'].map(model_files_dict) 

#delete extra columns
cleaned_otus_with_models.drop(['species_matched_model_list', 'genus_matched_model_list', 
                               'family_matched_model_list', 'model_columns_list'], axis=1, inplace=True)

#there is probably an easier way to do this, but I don't know it. 

#########################Part 2: merging models#################

#using a micom command:
#genus_model = micom.uti.join_models(["path_to_model_1.xml","path_to_model_2.xml","path_to_3.xml"], "cool_genus_model”).

#format should be otu_[number]_model = microm.uti.join_models(["model1.xml", "model2.xml",], "otu_[number]_model.xml")

#basically, for each row in otus_with_models, want to create a model from the values in model_columns_list
#and assign it to an OTU number. May be best to do this with a stupid dictionary. 

#for each row in otus_with_models [model_columns_list]:
#   genus_model = micom.uti.join_models([row], "cool_genus_model")
#   combo_model_list.append(genus_model)    

combined_model_list = []

#cleaned_otus_with_models2 = copy.deepcopy(cleaned_otus_with_models)
cleaned_otus_with_models3 = copy.deepcopy(cleaned_otus_with_models)


#is it a list vs string issue?
#cleaned_otus_with_models2['model_file_list_strings'] = ['", "'.join(map(str, l)) for l in cleaned_otus_with_models2['model_file_list']]
cleaned_otus_with_models3['model_file_list_strings'] = [', '.join(map(str, l)) for l in cleaned_otus_with_models3['model_file_list']]

#yes, it is a string issue. need to ahve files in string form and also in quotes.
#can't get the latter figured out. The above one creates a new column and separates the elements with ", ". But it doesn't
#put quotes in front of the first and last elements. 
#cleaned_otus_with_models2.update(cleaned_otus_with_models2[['model_file_list_strings']].applymap('"{}"'.format))    
#ok, tht worked

#cleaned_otus_with_models3['model_file_list_strings'].to_csv("00_model_strings.csv") #Just for checking

#this takes way too long to keep updating it every time there's a problem
#need to pull the file names from the stupid directory and compare them to model_file_list_strings

from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining") if isfile(join("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining", f))]

modelfiles = []
for elem in onlyfiles:
    if elem.endswith(".xml"):
        modelfiles.append(elem)

modelfiles_to_combine = []

for row in cleaned_otus_with_models3['model_file_list_strings']:
    rowlist = row.split(", ")
    for elem in rowlist:
        print(elem)
        if elem not in modelfiles_to_combine:
            modelfiles_to_combine.append(elem)

missing_models = np.setdiff1d(modelfiles_to_combine, modelfiles)
# yields the elements in `list_2` that are NOT in `list_1`
#only four models. Cool. 
#at this point I give up and I am editing the model file names
#notes: 
    #Clostridium_botulinum_A_str_Hall.xml becomes Clostridium_botulinum_A_str_Ha.xml
    #Clostridium_botulinum_E1_str_BoNT_E_Beluga.xml is in the model list as Clostridium_botulinum_E1_str_’BoNT_E_Beluga’.xml: pull out ‘ in .replace code
    #Clostridium_saccharoperbutylacetonicum_N1_4HMT.xml: this one is FIXED: used to be an underscore between the 4 and the H
    #Serratia_marcescens_subsp_marcescens_Db11.xml is old; now it is Serratia_marcescens_Db11.xml
    #any file that's _.xml I edited to be just .xml


#this will theoretically run join_models.    
for row in cleaned_otus_with_models3['model_file_list_strings']:
    #print(row)
    #print("      ")
    model_list = row.split(", ")
    #print(model_list)
    otu_model = micom.util.join_models(model_list, "otu_model")
    combined_model_list.append(otu_model)
#so, that took longer to time out than the others. Stalled out at Cap Capnocytophaga_ochracea_DSM_7271_.xml. I may just go in and hand edit that one.
#great news! My life is now incrimentally correcting the issues with specific file names in this code. 

#Well, that join_models code took, I shit you not, an hour. Hopefully it worked!
    
#now combined_model_list contains a model for each OTU in our dataset. Let's make a dictionary of them to dunk    
    
#make a dictionary where otus are matched with their lists of associated model files
model_files_dict = dict(zip(otu_list, combined_model_list)) #good thing lists are ordered, that's all I'll say. 

#so, now the trick is to iterate through the dictionary and save each model as an xml file. 
for key, value in model_files_dict.items():
    filename = "otu_"+str(key)+"model.xml"
    cobra.io.write_sbml_model(value, filename)  #I don't know if this will work, but maybe?
    
#start 7:05pm. End 8:30pm. Thank fuck the post-join models code worked. 


#dunk this as a new column in cleaned_otus_with_models
cleaned_otus_with_models['model_file_list'] = cleaned_otus_with_models['OTU_ID'].map(model_files_dict) 

        

#######testing grounds#################testing grounds#################testing grounds#################testing grounds##########
#######testing grounds#################testing grounds#################testing grounds#################testing grounds##########
#######testing grounds#################testing grounds#################testing grounds#################testing grounds##########

#should be able to pilot this with the MAMBO models
#Butyricimonas (genus, OTU 1135042) and Odoribacter(genus, OTU 196131) each have two MAMBO models assigned to them. 

#hello, I work            
otu_196131_model = micom.util.join_models(["Odoribacter_0.xml", "Odoribacter_1.xml"], "otu_196131_model.xml")
cobra.io.write_sbml_model(otu_196131_model, "otu_196131_model.xml")

#probably going to have to write some sort of complex function.
#first, need to add .xml extensions onto EVERY model name in every list in model_columns_list
#probably create a new column called model_files_list

#df['Discounted_Price'] = df.apply(lambda row: row.Cost - (row.Cost * 0.1), axis = 1) 


for row in otus_with_models['model_columns_list']:
    for elem in row:
        #print(elem)
        elem2 =elem + ".xml"
        print(elem2)

####################make a LIST of all unique models###############
#to get a list of all models, concatenate all models from all rows into one BIG list, then find unique values
longform_model_list = []

for row in otus_with_models['model_columns_list']:
    for elem in row:
        longform_model_list.append(elem)
        
unique_model_list = []

for elem in longform_model_list:
    print(elem)
    if elem not in unique_model_list:
        unique_model_list.append(elem)

unique_model_df = pd.DataFrame(unique_model_list)
unique_model_df.to_csv("04.30.20_unique_models_to_download.csv")        

#today's project: download ALL those fun models.

  
 
    

