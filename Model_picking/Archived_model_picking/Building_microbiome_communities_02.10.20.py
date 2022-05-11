#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 10:29:36 2020

@author: adamo010
"""
#in a classic turn of weirdness, had to install every fucking package
#very likely this is because I'm not running Spyder through Anaconda, which
#is because Spyder won't run through Anaconda anymore for some unknown reason

#anyway, to install gurobi:
#step 0: in terminal: pip install micom
#step 0.1: in terminal: install gurobi
#conda config --add channels http://conda.anaconda.org/gurobi
#conda install gurobi
#grbgetkey ed7e7426-3d57-11ea-a2f6-0a7c4f30bdbe
#(this is the academic license I signed up for on 01.22.20)
#note that I had to import it into the virtual environment to get it to run

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

#now that we have MICOM running and some real data, need to figure out how to build communities.
#I think the pipeline will be as follows:
    #import pandas series of phylogenetic breakdown from Burns OTU table and AGORA model list
    #add a new column in Burns OTU table where each cell is a(n empty for now) list
    #for each (start with genus, then go to family) in the Burns OTU table, find all models with a matching genus in
    #the AGORA model list, and append those models to the new list column in the Burns OTU table

#But right now Spyder isn't working so will probably have to trash my computer anyway. 
#fixed. Seriously, I hate virtual environments

agora = pd.read_csv("AGORA_model_phylogeny_table.csv")
otus = pd.read_csv("Burns_2015_OTU_phylogeny_table.csv")

#step 1: add a new column of empty lists to otus UPDATE DON'T ACTUALLY NEED THIS
#otus['model_list'] = np.empty((len(otus), 0)).tolist()
#for row in otus['model_list']:
    #print(type(row))

#step 2: add models to the model list YIKES

agora_genus_dict = agora.set_index('genus').to_dict()['model_ID']
#this creates a dictionary where the genus is the key and model_ID is the value
#NOPE this won't work, because every key has to be different and there are duplicate genuses.
agora_genus_dict2 = agora.set_index('model_ID').to_dict()['genus']
#THIS creates a dictionary where the model_ID is the value and genus is the key
#hopefully this will make the lookup easier.    
#I don't know if this is useful

#maybe I can pre-create genus lists of models. If genus is shared, append model_id to the genus list
unique_agora_genera = agora.genus.unique()

unique_agora_genera2 = pd.DataFrame(unique_agora_genera)    #import to pandas
unique_agora_genera2['model_list'] = np.empty((len(unique_agora_genera2), 0)).tolist()  #add an empty 'model list' column
unique_agora_genera2.reset_index()
unique_agora_genera2= unique_agora_genera2.rename(columns={0: "genus"})       #rename column

#what I want to do:
#for each row in unique_agora_genera2:
#scan through the genus column in agora
#if the value in the genus column in agora matches the value in the genus column of unique_agora_genera2,
#append the value of model_id in the agora dataframe to the list in the model_list column associated with that genus

models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]

for row in unique_agora_genera2['genus']:  
    genus_name_key = str(row)
    temp_df = agora[agora["genus"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    #now, I want to pull out the model names from these and append them to unique_agora_genera2       
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    models_dict[genus_name_key] = templist
    del templist    
    del genus_name_key

#EAT MY SHORTS PYTHON it finally worked. Precious models_dict
#now, do a better job of specifying different levels
genus_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_agora_genera2['genus']:  
    genus_name_key = str(row)
    temp_df = agora[agora["genus"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    #now, I want to pull out the model names from these and append them to unique_agora_genera2       
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    genus_models_dict[genus_name_key] = templist
    del templist    
    del genus_name_key
            
genus2_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_agora_genera2['genus']:  
    genus_name_key = str(row)
    temp_df = agora[agora["genus2"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    #now, I want to pull out the model names from these and append them to unique_agora_genera2       
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    genus2_models_dict[genus_name_key] = templist
    del templist    
    del genus_name_key
                
#interesting- looks like genus vs genus2 (AGORA-exact or renamed to match new species names, respectively) don't matter that much
    
#now for the horrid part: matching to the otu data in otus

model_matched_otus = copy.deepcopy(otus)    

model_matched_otus['genus_matched_model_list'] = model_matched_otus['genus'].map(genus_models_dict)        
        
#actually, no problem because I'm great.
#let's take a look in excel because I'm actually not that great
model_matched_otus.to_csv("Burns_2015_OTUs_matched_models.csv")

############################################################################################
#OK, it's a new day and I have a new plan
#need to do this iteratively:
#match models to taxa with species-level data
#pull those taxa out of the OTUs list
#pull out taxa WITH species-level information that LACK a model- may make models for these
#match models to (remaning) taxa with genus-level data (lists of matching models required)
#pull those taxa out of the OTUs list
#pull out taxa WITH genus-level information that LACK a model- may make models for these
#match models to (remaning) taxa with family-level data (lists of matching models required)
#pull those taxa out of the OTUs list
#pull out taxa WITH family level information that LACK a model- may make models for these
#see where we're at as far as coverage

#part 0: remove all obvious contaminants
#see OTU_table_filtering_pipeline_development.py






#################### JUNK
for row in unique_agora_genera2['genus']:
    #print(row)
    if str(row) in agora['genus']:
        print(row)



for row in otus['genus']:         #for each row in the column taxon_ID in the dataframe otus:
   print(row)
   if row =  
    
    
    model_list_dict.update({str(row),[]})

#testing     
    temp_df = agora[agora["genus"].str.contains('Aeromonas')]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    #now, I want to pull out the model names from these and append them to unique_agora_genera2       
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    models_dict['Aeromonas'] = templist
    del temp_df
    del templist    
        


for index, row in filtered_otus.iterrows():   
    col_list = colnames
    counter2 = 0
    if col in colnames:
        print(col)
        for col in col_list:
            if row >= 0.001:
                counter2 = counter2+1
            fitered_otus["num_above_cutoff"] = counter2    
    
#LowAb_otus = preprocessed_otus[preprocessed_otus]        
#preprocessed_otus.drop(LowAb_otus, inplace=True)   

#for col in preprocessed_otus.columns:
    #counter = 0
    #if 'Sample_' in str(col):
       # for row in col:
          #  if row >=0.001:
           #     counter= counter+1
    i#f counter >=18:
#












    
    