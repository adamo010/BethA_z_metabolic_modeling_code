#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 14:38:38 2020

@author: adamo010
"""

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


#all right, we have a list of OTUs with AGORA models (03.27.20_OTUs_with_models_V1.csv and 03.27.20_OTUs_with_models_V2.csv)
#we also have a list of OTUs without models (03.27.20_OTUs_without_models_V1.csv and 03.27.20_OTUs_without_models_V2.csv)
#in both of the without-model lists, there are 504 OTUs identified to at least the family level; 252 OTUs are at higher levels
#what other databases are available?

#The MAMBO models are available: https://github.com/danielriosgarza/MAMBO
#downloaded these: let's take a look: see MAMBO_models
#unlike stupid AGORA, can actually download all these models at once. 

#I am so lucky to have such resources at my disposal. MAMBO gives me all the models at once, but no taxonomy.
#AGORA gives me taxonomy, but I have to download all the models individually.
#what a time to be alive.

#print the first line of each file in the directory

with open("Dietzia_0.xml") as myfile:
    head =[next(myfile) for x in range(3)]
print(head)  

#Only some (probably few) of the models have their SPECIES listed, let alone a taxonomy. What the fuck. How is this database useable?

#just make a list of the file names- these are probably genera
from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir("/Users/adamo010/Documents/MICOM_CRC_FBA/MAMBO_models") if isfile(join("/Users/adamo010/Documents/MICOM_CRC_FBA/MAMBO_models", f))]

onlyfiles
model_genera = []
for modelname in onlyfiles:
    #print(type(modelname))
    #modelname.replace(".xml", "")
    #print(modelname)
    #modelname = modelname[:modelname.find('_')+1]
    modelname2 = modelname.split('_')[0]
    model_genera.append(modelname2)

model_genera2 = pd.Series(model_genera)     #convert list to pandas series
       
model_genera2.to_csv("MAMBO_model_genera.csv")    #save as csv: this is the genus?

#can we just match the genus-species name? Since getting family won't be trivial? Or at least just get genus species so I 
#can look up family by hand?

#would really love to pull model_id, if I can. third line of each file.
#google recommended putting this into a terminal window: head -n3 *.xml > MAMBO_model_ids.txt
#it mostly worked. now I have a text file called MAMBO_model_ids.txt that contains the first three lines of each model file,
#plus the model name. Which is a step forward, I guess.

MAMBO_model_list_key = []
MAMBO_model_list_value = []

with open("MAMBO_model_ids.txt") as myfile:
    for line in myfile:
        temp = str(line)    #have to do this b/c otherwise each element's type will be 'type' in the ensuing list. For some reason.
        if "==>" in temp:
            MAMBO_model_list_key.append(str(temp))
        if "model id" in temp:
            #print(line)
            MAMBO_model_list_value.append(str(temp))

#now I have two lists: one with the model name and one with the model_id (which is kind of like species? in some cases?)
#need to clean these up.

#this code splits up the model_list_key into a list of space-separated elements, pulls out the actual model file name,
#and sticks it in a new list MAMBO_model_list_key2            
MAMBO_model_list_key2 = []
for elem in MAMBO_model_list_key:
    elem2 = elem.split()
    #print(elem2)
    for subelem in elem2:
        if '.xml' in subelem:
            MAMBO_model_list_key2.append(subelem)

#try something similar for list_value
MAMBO_model_list_value2 = []
for elem in MAMBO_model_list_value:
    elem2 = elem.split()
    #print(elem2)
    for subelem in elem2:
        if 'id=' in subelem:
            subelem= subelem.lstrip("id=")
            subelem= subelem.strip('"')
            MAMBO_model_list_value2.append(subelem)

MAMBO_dictionary = dict(zip(MAMBO_model_list_key2, MAMBO_model_list_value2))             
#okay, so for some of these it LOOKS like the model and the filename aren't matching up.
#e.g. for Acetomicrobium_1.xml, the model_id is Anaerobaculum hydrogeniformans.
#this is fine. The .xml name is the new genus. or an updated genus. something like that. Real helpful, right?

#all right. Next, we want to get this dictionary into a pandas series. Then, we want to extract the species names from the values.
#then, I guess we hand-look up all the phylogenies of these fuckers. can probably collapse at the species level. Hopefully.

#values that start with A_ or kb_ are not helpful re: species names

mambo_model_table= pd.DataFrame(list(MAMBO_dictionary.items()),columns = ['model_file_name','model_id'])

cursed_elems = ["kb", "A"]

#for elem in mambo_model_table["model_id"]:
  #  namelists = elem.split("_")
    #if namelists[0] not in cursed_elems:
      #  print(namelists)
      #  genus = namelists[1]
       # species = namelists[2]

mambo_model_table['model_id_list']= mambo_model_table['model_id'].str.split("_")    

genus_list = []
species_list = []   #I hope no one ever has to use this code since it's so hacky.

for row in mambo_model_table['model_id_list']:
    genus = str(row[0])
    species = str(row[1])
    genus_list.append(genus)
    species_list.append(species)

mambo_model_table['Genus'] = genus_list
mambo_model_table['Species'] = species_list
    
mambo_model_table["Genus_species"] = mambo_model_table['Genus'] + " " + mambo_model_table['Species']

mambo_model_table.to_csv("03.20.20_MAMBO_model_spp_for_id.csv")

#now for the boring work.
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
