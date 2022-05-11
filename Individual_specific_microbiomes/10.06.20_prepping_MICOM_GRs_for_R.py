#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 19:59:20 2020

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

results_file_list =[]
for file in os.listdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/"):
    src= file
    if fnmatch.fnmatch(file, "*_fba_results_plus_taxonomies.csv"):
        results_file_list.append(str(file))
del file
del src        

sample_name_list = []
for elem in results_file_list:
    elem2 = elem.rstrip("_fba_results_plus_taxonomies.csv")
    sample_name_list.append(elem2)
del elem
del elem2       

results_file_dict = dict(zip(sample_name_list, results_file_list)) 

data_table_list = []
data_table_names = []

for key, value in results_file_dict.items():
    file = pd.read_csv(value)
    data_table_list.append(file)
    data_table_names.append(key)
del key
del value    

data_file_dict = dict(zip(data_table_names, data_table_list)) 

#now have all tables imported. Each table has the following headings:
#OTU_ID, abundance, growth_rate, reactions, metabolites, taxonomy: all are numerical except for taxonomy

#look at fusobacterium, faecalibacterium, Bfrag specifically?
#look at top 10% of growers in each community?

#fusobacterium OTU_ID: 1096339; e. coli 656881;  Bfrag: 128205; peptostreptococcus: 3878210
#Protective: Faecalibacterium: 147702; lactobacillus: 107784

#I would like to pull specific OTUs from each of the 88 sample dataframes (abundance and growth_rate), 
#and add a new column called "sample_ID"

Fuso_full_df = pd.DataFrame(columns=["OTU_ID", "abundance", "growth_rate", "Sample_ID"])
Ecoli_full_df = pd.DataFrame(columns=["OTU_ID", "abundance", "growth_rate", "Sample_ID"])
Faecal_full_df = pd.DataFrame(columns=["OTU_ID", "abundance", "growth_rate", "Sample_ID"])
Lacto_full_df = pd.DataFrame(columns=["OTU_ID", "abundance", "growth_rate", "Sample_ID"])

for key, value in data_file_dict.items():
    value['OTU_ID'] = value['OTU_ID'].astype(str) 
    value.set_index("OTU_ID")
    Fuso_df= value[value['OTU_ID'] == "1096339"]
    Fuso_df2= copy.deepcopy(Fuso_df)
    Fuso_df2['Sample_ID'] = str(key)
    Fuso_full_df= Fuso_full_df.append(Fuso_df2)
    Ecoli_df= value[value['OTU_ID'] == "656881"]
    Ecoli_df2= copy.deepcopy(Ecoli_df)
    Ecoli_df2['Sample_ID'] = str(key)
    Ecoli_full_df= Ecoli_full_df.append(Ecoli_df2)
    Faecal_df= value[value['OTU_ID'] == "147702"]
    Faecal_df2= copy.deepcopy(Faecal_df)
    Faecal_df2['Sample_ID'] = str(key)
    Faecal_full_df= Faecal_full_df.append(Faecal_df2)
    Lacto_df= value[value['OTU_ID'] == "107784"]
    Lacto_df2= copy.deepcopy(Lacto_df)
    Lacto_df2['Sample_ID'] = str(key)
    Lacto_full_df= Lacto_full_df.append(Lacto_df2)
    
#oh cool, every GR is the same for every community.
    
########    
#pull top 10 most abundant otus from each dataframe
growth_rate_df = pd.read_csv("10.07.20_OTU_ids_only.csv")  
growth_rate_df['OTU_ID']=growth_rate_df['OTU_ID'].astype(str) 
    
for key, value in data_file_dict.items():
    top10 = value.nlargest(15,'growth_rate') 
    top10_filtered = top10[['OTU_ID', 'growth_rate']]
    top10_filtered2 = copy.deepcopy(top10_filtered)
    top10_filtered2['OTU_ID']=top10_filtered2['OTU_ID'].astype(str) 
    top10_filtered2.rename(columns = {'growth_rate':key}, inplace = True) 
    #print(top10_filtered)
    growth_rate_df = pd.merge(growth_rate_df, top10_filtered2, how="left", left_on="OTU_ID", right_on="OTU_ID")
    
growth_rate_df.to_csv("10.07.20_top_15_GRs_from_MICOM.csv")    






