#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 11:22:17 2021

@author: adamo010
"""

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
import pickle
import sys #for specifying input files externally


#the goal of this code is to prepare the growth rate outputs from MICOM on MSI for procesisng in R.
#will probably be some back and forth because, let's face it, R sucks and I often think I can code stuff there (esp wrt dataframe
#manipulation) and it goes horribly wrong and I end up back here. With ol' python.

#anyway.
#based on 01.22.21_EUstd_MM_comm_and_indiv_GR_processing.py and 08.31.21_Hale2018_comm_GR_combining_files_and_adding_metadata.py

#let's start by bringing in some files
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Common_growth_medium_runs/June2021_rerunning_MICOM_on_MSI_with_fluxes_on/Sept2021_V30_output")

comm_GR_file_list = glob.glob('*_V30_Burns_fluxes_on_comm_GRs.csv')

comm_GR_file_names = []

for elem in comm_GR_file_list:
    sample_name = re.sub("_V30_Burns_fluxes_on_comm_GRs.csv", "", elem)  
    comm_GR_file_names.append(sample_name)   
       
del elem, sample_name    

#START HERE
comm_GRs =  {}
for elem in comm_GR_file_list:
    output_file = pd.read_csv(elem)
    sample_name = re.sub("_V30_Burns_fluxes_on_comm_GRs.csv", "", elem) #delete extra crap on string.
    #print(sample_name)
    for column in output_file:
        costring1 = re.sub("<CommunitySolution ", "", column) #delete first part of string from community solution
        print(costring1)
        costring2 = costring1[:-19] #delete last 19 characters of string. These all differ, so can't match text/symbols
        print(costring2)
        comm_GRs.update({sample_name: costring2})
del elem, column, costring1, costring2, output_file, sample_name

#nice. Save our output file

os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/") #analysis directory
with open('Burns2015_V30MSIrun_community_GRs_combined.csv', 'w') as f:
    for key in comm_GRs.keys():
        f.write("%s,%s\n"%(key,comm_GRs[key]))
del f, key, comm_GR_file_list, comm_GRs, comm_GR_file_names        
        
#NOWWWWW we want to make something that we can use with the R code I wrote, so I don't have to re-write the R code.
#Step 1: importing the dataframes. 
TD_point_two = pd.read_csv("Burns2015_V30MSIrun_community_GRs_combined.csv", names= ["Sample_ID", "comm_gr"])
#TD isn't point two. It's point one. But who cares. 

metadata = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Burns_2015_sample_metadata.csv")
#wtf, why are therea bunch of extra rows here?
metadata.dropna(subset = ["Patient_Blind_ID"], inplace=True)

comm_GRs_with_metadata = pd.merge(left=TD_point_two, right=metadata, how='left', left_on='Sample_ID', right_on='SampleID')
comm_GRs_with_metadata.drop(columns=['SampleID'], inplace=True)

#for another analysis, I would like to convert my data (non-summarized) from wide to long format.
#Python, of course, has a function for this.
#but first, set column names as strings
comm_GRs_with_metadata.columns = comm_GRs_with_metadata.columns.astype(str)
#reset comm_gr as float.
comm_GRs_with_metadata["comm_gr"] = comm_GRs_with_metadata["comm_gr"].astype(float)
#add rankings
comm_GRs_with_metadata["rank"] = comm_GRs_with_metadata.groupby(["Description"])["comm_gr"].rank(method='first', ascending=False)
#IMPORTANT NOTE: rank method=first is a way to break ties. This method specifies that the first incidence of a value gets the higher rank.
#method=dense will... possibly drop samples? Unclear. When I used dense, it dropped a bunch of data for some inexplicable reason

comm_GRs_with_metadata.to_csv("09.09.21_Burns2015_V30MSIrun_community_GRs_combined_plus_ranks.csv")





