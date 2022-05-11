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


#all right, we have a list of OTUs with AGORA models (03.18.20_OTUs_with_models_V1.csv and 03.18.20_OTUs_with_models_V2.csv)
#we also have a list of OTUs without models (03.18.20_OTUs_without_models_V1.csv and 03.18.20_OTUs_without_models_V2.csv)
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
       
model_genera2.to_csv("MAMBO_model_genera.csv")    #save as csv

#model picking code is missing a bunch of fucking models for some inexplicable reason. Have to go revisit that.     


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
