#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 16:48:43 2021

@author: adamo010
"""
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
import pickle
import sys #for specifying input files externally


#the goal of this code is to prepare the growth rate outputs from MICOM on MSI for procesisng in R.
#will probably be some back and forth because, let's face it, R sucks and I often think I can code stuff there (esp wrt dataframe
#manipulation) and it goes horribly wrong and I end up back here. With ol' python.

#anyway.
#based on 01.22.21_EUstd_MM_comm_and_indiv_GR_processing.py
#and 08.31.21_Hale2018_comm_GR_combining_files_and_adding_metadata.py

#let's start by bringing in some files
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Muehlbauer2021/Muehlbauer2021_output_from_MSI/")

comm_GR_file_list = glob.glob('*_V01_Muehlbauer_comm_GRs.csv')

comm_GR_file_names = []

for elem in comm_GR_file_list:
    sample_name = re.sub("_V01_Muehlbauer_comm_GRs.csv", "", elem)  
    comm_GR_file_names.append(sample_name)        
del elem, sample_name    

#START HERE
comm_GRs =  {}
for elem in comm_GR_file_list:
    output_file = pd.read_csv(elem)
    sample_name = re.sub("_V01_Muehlbauer_comm_GRs.csv", "", elem) #delete extra crap on string.
    for column in output_file:
        costring1 = re.sub("<CommunitySolution ", "", column) #delete first part of string from community solution
        #print(costring1)
        costring2 = costring1[:-19] #delete last 19 characters of string. These all differ, so can't match text/symbols
        #print(costring2)
        comm_GRs.update({sample_name: costring2})
del elem, column, costring1, costring2, output_file, sample_name

#nice. Save our output file

os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Muehlbauer2021/")

with open('Muehlbauer2021_community_GRs_combined.csv', 'w') as f:
    for key in comm_GRs.keys():
        f.write("%s,%s\n"%(key,comm_GRs[key]))
del f, key, comm_GR_file_list, comm_GRs, comm_GR_file_names        
        
#NOWWWWW we want to make something that we can use with the R code I wrote, so I don't have to re-write the R code.
#Step 1: importing the dataframes. 
TD_point_two = pd.read_csv("Muehlbauer2021_community_GRs_combined.csv", names= ["Sample_ID", "comm_gr"])

metadata = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/ASV_tables_Muehlbauer2021/Muehlbauer2021_metadata_16s_only.csv")

comm_GRs_with_metadata = pd.merge(left=TD_point_two, right=metadata, how='left', left_on='Sample_ID', right_on='Sample_name')
comm_GRs_with_metadata.drop(columns=['Sample_name'], inplace=True)

#for another analysis, I would like to convert my data (non-summarized) from wide to long format.
#Python, of course, has a function for this.
#but first, set column names as strings
comm_GRs_with_metadata.columns = comm_GRs_with_metadata.columns.astype(str)
#I'm not sure if we need this right now, but hang onto it just in case. I only ran one tradeoff for Hale data, so no need to
#convert to long. 

#also, set column names to something with a call-able string; in the pd.wide_to_long code, a common character is required
#so Python knows which columns it needs to shuffle together.

#comm_GRs_with_metadata.rename(columns={"0.5": "0.5_tdof", "0.25": "0.25_tdof", "0.1": "0.1_tdof"}, inplace = True)

#EDIT ME!!!!!!! I'm stacking extra columns weirdly.
#comm_GRs_with_metadata.drop(columns=['Site', "MSI_status", "Stage"], inplace=True)
#longform_comm_GRs = pd.wide_to_long(comm_GRs_with_metadata, stubnames ="_tdof", i= ["Sample_ID", "Description", "Patient_Blind_ID"], 
                                    #j="tradeoff", sep = "_")

#longform_comm_GRs = comm_GRs_with_metadata.melt(id_vars=["Sample_ID", "Description", "Patient_Blind_ID"], 
                              #var_name="tradeoff", 
                              #value_name="_tdof")

#longform_comm_GRs= longform_comm_GRs.rename(columns={"_tdof":"Growth_rate"})
#longform_comm_GRs['tradeoff'] = longform_comm_GRs['tradeoff'].map(lambda x: x.rstrip('_tdof'))

#longform_comm_GRs.to_csv("01.25.21_EUstd_diet_longform_comm_GRs_by_tradeoff_plus_metadata.csv")

#adding rankings
#also not sure if I need this for Hale data, but can always drop it later
comm_GRs_with_metadata["comm_gr"] = comm_GRs_with_metadata["comm_gr"].astype(float)

comm_GRs_with_metadata["rank"] = comm_GRs_with_metadata.groupby(["experimental_treatment"])["comm_gr"].rank(method='first', ascending=False)
#IMPORTANT NOTE: rank method=first is a way to break ties. This method specifies that the first incidence of a value gets the higher rank.
#method=dense will... possibly drop samples? Unclear. When I used dense, it dropped a bunch of data for some inexplicable reason

comm_GRs_with_metadata.to_csv("12.29.21_Muehlbauer2021_community_GRs_combined_plus_ranks.csv")

#OKAY!!!!! Now we do the individual species growth rate stuff, but assign taxonomy to it first. 
#indivspp_GR_file_list = glob.glob('*_indiv_spp_GRs_EUstandard_MM_Jan2021.csv')
#indivspp_GR_file_names =[]
#for elem in indivspp_GR_file_list:
    #sample_name = re.sub("_indiv_spp_GRs_EUstandard_MM_Jan2021.csv", "", elem)  
    #indivspp_GR_file_names.append(sample_name)           
        
        
###################for other tradeoffs        
#indivspp_GR_file_list = glob.glob('*_indiv_spp_GRs_EUstandard_MM_Jan2021_point_two_five.csv')
#indivspp_GR_file_names =[]
#for elem in indivspp_GR_file_list:
    #sample_name = re.sub("_indiv_spp_GRs_EUstandard_MM_Jan2021_point_two_five.csv", "", elem)  
    #indivspp_GR_file_names.append(sample_name)       
    
#indivspp_GR_file_list = glob.glob('*_indiv_spp_GRs_EUstandard_MM_Jan2021_point_one.csv')
#indivspp_GR_file_names =[]
#for elem in indivspp_GR_file_list:
    #sample_name = re.sub("_indiv_spp_GRs_EUstandard_MM_Jan2021_point_one.csv", "", elem)  
    #indivspp_GR_file_names.append(sample_name)       







      