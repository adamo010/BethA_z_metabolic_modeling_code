#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 13:42:35 2020

@author: adamo010
"""

#Borrowed from 08.10.20_building_patient_specific_MICOM_comms.py
#goal is to count the number of models/OTUs within each patient MICOM community

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


os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA")

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

#strip the table to something useful.

otu_table_slimmed = otu_table_relabund.drop(columns=["taxonomy", "model_file_name"])

#create a new row that contains a count of all nonzero values in the column
for column in otu_table_slimmed:
    count = np.count_nonzero(column)

count_row = np.count_nonzero(otu_table_slimmed, axis=1)
count_row2= np.count_nonzero(otu_table_slimmed, axis=0) #this is the one I want.

#now, to get a CSV file with a sample_ID column and a count column
sample_IDs = otu_table_slimmed.columns.tolist()

count_dict = dict(zip(sample_IDs, count_row2))

with open('12.03.20_counts_num_models_per_MICOM_community.csv', 'w') as csv_file:  
    writer = csv.writer(csv_file)
    for key, value in count_dict.items():
       writer.writerow([key, value])
