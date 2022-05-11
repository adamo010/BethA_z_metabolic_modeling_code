#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:45:44 2020

@author: adamo010
"""

#updates from 09.28.20: shit won't run on MSI and it is MADDENING

#Have made community pickle files. Time to run some GROWTH RATES. Ideally, this would be on MSI, but since gurobi is apparently 
#unworkable on MSI without significant troubleshooting, let's try to run a few of these locally.
#use the directory /Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes
#copy most of this from 08.27.20_patient_specific_community_GRs.py

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
from micom import load_pickle
from micom.media import minimal_medium

#It's worth noting that I'd normally for-loop the hell out of this. But pickle files take so long to input that I'm not going to do that.
#What I might do is dump this into a function where the input is a pickle file. Then I can run each pickle file individually.

#step 0: make a list of all the pickle files

pickle_file_list = []

for file in os.listdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/"):
    src= file
    if fnmatch.fnmatch(file, "*_community.pickle"):
        pickle_file_list.append(str(file))
del file
del src        

#pickle_file_list now contains all the pickle files
sample_name_list = []
for elem in pickle_file_list:
    elem2 = elem.rstrip("_community.pickle")
    sample_name_list.append(elem2)
del elem
del elem2        

pickle_file_dict = dict(zip(sample_name_list, pickle_file_list))

#create four dictionaries to chunk this out so it can run locally

pickle_file_subdict_1 = {}
pickle_file_subdict_2 = {}
pickle_file_subdict_3 = {}
pickle_file_subdict_4 = {}

for k, v in pickle_file_dict.items():
    if len(pickle_file_subdict_1.keys()) < 22:
        pickle_file_subdict_1.update({k: v})
    elif len(pickle_file_subdict_2.keys()) < 22:
        pickle_file_subdict_2.update({k: v})
    elif len(pickle_file_subdict_3.keys()) < 22:
        pickle_file_subdict_3.update({k: v})
    elif len(pickle_file_subdict_4.keys()) < 22:
        pickle_file_subdict_4.update({k: v})    

community_output_list1 = []
community_name_list1 = []
community_output_list2 = []
community_name_list2 = []
community_output_list3 = []
community_name_list3 = []
community_output_list4 = []
community_name_list4 = []

        
for key, value in pickle_file_subdict_1.items():
    sample_comm = load_pickle(value) #takes a long time
    sample_sol = sample_comm.optimize() #short-ish (10 seconds)
    community_output_list1.append(sample_sol)
    community_name_list1.append(key)
    del sample_comm
    del sample_sol
os.system('say "Another one bites the dust"')

community_dict1 = dict(zip(community_name_list1, community_output_list1)) 
for key, value in community_dict1.items():
    print(key, value)

f = open("10.06.20_community_dict1.csv","w")
f.write( str(community_dict1) )
f.close()    

with open('10.06.20_community_dict1.csv', 'w') as csv_file:  
    writer = csv.writer(csv_file)
    for key, value in community_dict1.items():
       writer.writerow([key, value])

#now, I would like to save fluxes:
for key, value in community_dict1.items():
    print(value.fluxes)

#step 2
for key, value in pickle_file_subdict_2.items():
    sample_comm = load_pickle(value) #takes a long time
    sample_sol = sample_comm.optimize() #short-ish (10 seconds)
    community_output_list2.append(sample_sol)
    community_name_list2.append(key)
    del sample_comm
    del sample_sol
os.system('say "Another one bites the dust"')

community_dict2 = dict(zip(community_name_list2, community_output_list2)) 
for key, value in community_dict2.items():
    print(key, value)

with open('10.06.20_community_dict2.csv', 'w') as csv_file:  
    writer = csv.writer(csv_file)
    for key, value in community_dict2.items():
       writer.writerow([key, value])

#step 3
for key, value in pickle_file_subdict_3.items():
    sample_comm = load_pickle(value) #takes a long time
    sample_sol = sample_comm.optimize() #short-ish (10 seconds)
    community_output_list3.append(sample_sol)
    community_name_list3.append(key)
    del sample_comm
    del sample_sol
os.system('say "Another one bites the dust"')

community_dict3 = dict(zip(community_name_list3, community_output_list3)) 
for key, value in community_dict3.items():
    print(key, value)

with open('10.06.20_community_dict3.csv', 'w') as csv_file:  
    writer = csv.writer(csv_file)
    for key, value in community_dict3.items():
       writer.writerow([key, value])
        
 #step 4       
for key, value in pickle_file_subdict_4.items():
    sample_comm = load_pickle(value) #takes a long time
    sample_sol = sample_comm.optimize() #short-ish (10 seconds)
    community_output_list4.append(sample_sol)
    community_name_list4.append(key)
    del sample_comm
    del sample_sol
os.system('say "Another one bites the dust"')

community_dict4 = dict(zip(community_name_list4, community_output_list4)) 
for key, value in community_dict4.items():
    print(key, value)

with open('10.06.20_community_dict4.csv', 'w') as csv_file:  
    writer = csv.writer(csv_file)
    for key, value in community_dict4.items():
       writer.writerow([key, value])
        
########################update from 10.08.20###############
#Any given species has the same growth rate across all communities (although each species has a different growth rate)
#Ran suggested running each        
        
        
        