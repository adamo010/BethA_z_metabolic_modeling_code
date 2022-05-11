#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 13:01:25 2020

@author: adamo010
"""

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

#now I have a list of files, a list of sample names, and a dictionary combining both.

#step 1: write the function not in function form (yet)
from micom import load_pickle

sample_07_comm = load_pickle("Sample_07_community.pickle") #takes a long time
sample_07_sol = sample_07_comm.optimize() #short-ish (10 seconds)

print(sample_07_comm.objective.expression)
sample_07_comm.optimize()

GR = str(print(sample_07_sol))
growth_rate_subdict_1.update({key, GR})


#########the in between###########
from micom.media import minimal_medium
rates = sample_07_sol.members.growth_rate.drop("medium")
med = minimal_medium(sample_07_comm, sample_07_sol.growth_rate, min_growth=rates, exports=True) #oh this takes a long time too, huh?

#########the in between###########

output = open("_prefba_output_summary.csv", 'w')
results = sample_07_sol.members #this is a pandas dataframe
pd.set_option('display.max_rows', None)
output.write(str(results))
output.close()
del output
del results

with open("_fba_output_summary.csv", "w") as csvFile:        
    writer = csv.writer(csvFile)
    with open("_prefba_output_summary.csv") as input:
        fluxreader = csv.reader(input) 
        print(fluxreader)
        fluxreader2 = filter(None, fluxreader) #get rid of empty rows
        for row in fluxreader2:
            row = row[0].split()
            firstel = row[0]
            if len(row) < 5:
                row = ['compartments', 'abundance', 'growth_rate', 'reactions', 'metabolites']
            writer.writerow(row)                    
df = pd.read_csv('_fba_output_summary.csv', skiprows=1)
df.to_csv('_fba_output_summary.csv', index=False)
os.remove("_prefba_output_summary.csv")      
for file in os.listdir():                       
    src=file
    if fnmatch.fnmatch(file, "_fba_output_summary.csv"):
        dst = "Sample_07"+file
        os.rename(src,dst)
del csvFile
del df
del dst
del file
del firstel
del fluxreader
del fluxreader2
del src
del writer
del input
del row
      

#Step 2: getting taxonomies in
#okay, where is my dictionary of otus (keys) and taxonomies (values)?
preprocessed_otus = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal/Burns_absolute_OTU_table_phylogeny_corrected.csv")
preprocessed_otus.rename(columns={"#OTU_ID":"OTU_ID"}, inplace = True)
#only keep first and last column
otu_taxonomy_lookup_table = copy.deepcopy(preprocessed_otus)
otu_taxonomy_lookup_table = otu_taxonomy_lookup_table[['OTU_ID','taxonomy']]
#need to convert this to string for downstream stuff
otu_taxonomy_lookup_table["OTU_ID"] = otu_taxonomy_lookup_table["OTU_ID"].astype(str)

#step 3: make taxonomy-adapted CSV table
sample_07_pretax = pd.read_csv("Sample_07_fba_output_summary.csv")
sample_07_pretax.columns=["OTU_ID", "abundance", "growth_rate", "reactions", "metabolites"]

sample_07_output = pd.merge(sample_07_pretax, otu_taxonomy_lookup_table, on="OTU_ID", how="left")

sample_07_output.to_csv("_fba_results_plus_taxonomies.csv", index=False)
for file in os.listdir():
    src=file
    if fnmatch.fnmatch(file, "_fba_results_plus_taxonomies.csv"):
        dst = str("Sample_07")+file
        os.rename(src,dst)
del dst
del file
del src

#great. Sample 07 is done. Only 87 more to go. 

#smash this all into a function...

def growth_rates_from_FBA(file_name, sample_name):
    sample_comm = load_pickle(file_name) #takes a long time
    sample_sol = sample_comm.optimize() #short-ish (10 seconds)
    output = open("_prefba_output_summary.csv", 'w')
    results = sample_sol.members #this is a pandas dataframe
    pd.set_option('display.max_rows', None)
    output.write(str(results))
    output.close()
    del output
    del results
    with open("_fba_output_summary.csv", "w") as csvFile:        
        writer = csv.writer(csvFile)
        with open("_prefba_output_summary.csv") as input:
            fluxreader = csv.reader(input) 
            print(fluxreader)
            fluxreader2 = filter(None, fluxreader) #get rid of empty rows
            for row in fluxreader2:
                row = row[0].split()
                firstel = row[0]
                if len(row) < 5:
                    row = ['compartments', 'abundance', 'growth_rate', 'reactions', 'metabolites']
                writer.writerow(row)                    
    df = pd.read_csv('_fba_output_summary.csv', skiprows=1)
    df.to_csv('_fba_output_summary.csv', index=False)
    os.remove("_prefba_output_summary.csv")      
    for file in os.listdir():                       
        src=file
        if fnmatch.fnmatch(file, "_fba_output_summary.csv"):
            dst = str(sample_name)+file
            os.rename(src,dst)
    del csvFile
    del df
    del dst
    del file
    del firstel
    del fluxreader
    del fluxreader2
    del src
    del writer
    del input        
    return

def add_taxonomies_to_GRs(file_name, sample_name):
    sample_pretax = pd.read_csv(file_name)
    sample_pretax.columns=["OTU_ID", "abundance", "growth_rate", "reactions", "metabolites"]
    sample_output = pd.merge(sample_pretax, otu_taxonomy_lookup_table, on="OTU_ID", how="left")
    sample_output.to_csv("_fba_results_plus_taxonomies.csv", index=False)
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_fba_results_plus_taxonomies.csv"):
            dst = str(sample_name)+file
            os.rename(src,dst)
    del dst
    del file
    del src
    os.system('say "I have completed my assigned task your ladyship"')
    return

#try it with Sample 39 (next on the list)
growth_rates_from_FBA("Sample_39_community.pickle", "Sample_39")
add_taxonomies_to_GRs("Sample_39_fba_output_summary.csv", "Sample_39")
    
#works. Moving on...
growth_rates_from_FBA("Sample_84_community.pickle", "Sample_84")
add_taxonomies_to_GRs("Sample_84_fba_output_summary.csv", "Sample_84")

growth_rates_from_FBA("Sample_62_community.pickle", "Sample_62")
add_taxonomies_to_GRs("Sample_62_fba_output_summary.csv", "Sample_62")

growth_rates_from_FBA("Sample_43_community.pickle", "Sample_43")
add_taxonomies_to_GRs("Sample_43_fba_output_summary.csv", "Sample_43")

growth_rates_from_FBA("Sample_18_community.pickle", "Sample_18")
add_taxonomies_to_GRs("Sample_18_fba_output_summary.csv", "Sample_18")

growth_rates_from_FBA("Sample_26_community.pickle", "Sample_26")
add_taxonomies_to_GRs("Sample_26_fba_output_summary.csv", "Sample_26")

growth_rates_from_FBA("Sample_57_community.pickle", "Sample_57")
add_taxonomies_to_GRs("Sample_57_fba_output_summary.csv", "Sample_57")

growth_rates_from_FBA("Sample_69_community.pickle", "Sample_69")
add_taxonomies_to_GRs("Sample_69_fba_output_summary.csv", "Sample_69")

growth_rates_from_FBA("Sample_32_community.pickle", "Sample_32")
add_taxonomies_to_GRs("Sample_32_fba_output_summary.csv", "Sample_32")

growth_rates_from_FBA("Sample_13_community.pickle", "Sample_13")
add_taxonomies_to_GRs("Sample_13_fba_output_summary.csv", "Sample_13")

growth_rates_from_FBA("Sample_48_community.pickle", "Sample_48")
add_taxonomies_to_GRs("Sample_48_fba_output_summary.csv", "Sample_48")

growth_rates_from_FBA("Sample_76_community.pickle", "Sample_76")
add_taxonomies_to_GRs("Sample_76_fba_output_summary.csv", "Sample_76")

growth_rates_from_FBA("Sample_90_community.pickle", "Sample_90")
add_taxonomies_to_GRs("Sample_90_fba_output_summary.csv", "Sample_90")

growth_rates_from_FBA("Sample_29_community.pickle", "Sample_29")
add_taxonomies_to_GRs("Sample_29_fba_output_summary.csv", "Sample_29")

growth_rates_from_FBA("Sample_17_community.pickle", "Sample_17")
add_taxonomies_to_GRs("Sample_17_fba_output_summary.csv", "Sample_17")

growth_rates_from_FBA("Sample_72_community.pickle", "Sample_72")
add_taxonomies_to_GRs("Sample_72_fba_output_summary.csv", "Sample_72")

growth_rates_from_FBA("Sample_94_community.pickle", "Sample_94")
add_taxonomies_to_GRs("Sample_94_fba_output_summary.csv", "Sample_94")

growth_rates_from_FBA("Sample_53_community.pickle", "Sample_53")
add_taxonomies_to_GRs("Sample_53_fba_output_summary.csv", "Sample_53")

growth_rates_from_FBA("Sample_36_community.pickle", "Sample_36")
add_taxonomies_to_GRs("Sample_36_fba_output_summary.csv", "Sample_36")

growth_rates_from_FBA("Sample_08_community.pickle", "Sample_08")
add_taxonomies_to_GRs("Sample_08_fba_output_summary.csv", "Sample_08")

growth_rates_from_FBA("Sample_79_community.pickle", "Sample_79")
add_taxonomies_to_GRs("Sample_79_fba_output_summary.csv", "Sample_79")

growth_rates_from_FBA("Sample_47_community.pickle", "Sample_47")
add_taxonomies_to_GRs("Sample_47_fba_output_summary.csv", "Sample_47")

growth_rates_from_FBA("Sample_22_community.pickle", "Sample_22")
add_taxonomies_to_GRs("Sample_22_fba_output_summary.csv", "Sample_22")

growth_rates_from_FBA("Sample_80_community.pickle", "Sample_80")
add_taxonomies_to_GRs("Sample_80_fba_output_summary.csv", "Sample_80")

growth_rates_from_FBA("Sample_66_community.pickle", "Sample_66")
add_taxonomies_to_GRs("Sample_66_fba_output_summary.csv", "Sample_66")

growth_rates_from_FBA("Sample_58_community.pickle", "Sample_58")
add_taxonomies_to_GRs("Sample_58_fba_output_summary.csv", "Sample_58")

growth_rates_from_FBA("Sample_19_community.pickle", "Sample_19")
add_taxonomies_to_GRs("Sample_19_fba_output_summary.csv", "Sample_19")

growth_rates_from_FBA("Sample_27_community.pickle", "Sample_27")
add_taxonomies_to_GRs("Sample_27_fba_output_summary.csv", "Sample_27")

growth_rates_from_FBA("Sample_42_community.pickle", "Sample_42")
add_taxonomies_to_GRs("Sample_42_fba_output_summary.csv", "Sample_42")

growth_rates_from_FBA("Sample_63_community.pickle", "Sample_63")
add_taxonomies_to_GRs("Sample_63_fba_output_summary.csv", "Sample_63")

growth_rates_from_FBA("Sample_85_community.pickle", "Sample_85")
add_taxonomies_to_GRs("Sample_85_fba_output_summary.csv", "Sample_85")

growth_rates_from_FBA("Sample_06_community.pickle", "Sample_06")
add_taxonomies_to_GRs("Sample_06_fba_output_summary.csv", "Sample_06")

growth_rates_from_FBA("Sample_38_community.pickle", "Sample_38")
add_taxonomies_to_GRs("Sample_38_fba_output_summary.csv", "Sample_38")

growth_rates_from_FBA("Sample_91_community.pickle", "Sample_91")
add_taxonomies_to_GRs("Sample_91_fba_output_summary.csv", "Sample_91")

growth_rates_from_FBA("Sample_49_community.pickle", "Sample_49")
add_taxonomies_to_GRs("Sample_49_fba_output_summary.csv", "Sample_49")

growth_rates_from_FBA("Sample_77_community.pickle", "Sample_77")
add_taxonomies_to_GRs("Sample_77_fba_output_summary.csv", "Sample_77")

growth_rates_from_FBA("Sample_12_community.pickle", "Sample_12")
add_taxonomies_to_GRs("Sample_12_fba_output_summary.csv", "Sample_12")

growth_rates_from_FBA("Sample_33_community.pickle", "Sample_33")
add_taxonomies_to_GRs("Sample_33_fba_output_summary.csv", "Sample_33")

growth_rates_from_FBA("Sample_56_community.pickle", "Sample_56")
add_taxonomies_to_GRs("Sample_56_fba_output_summary.csv", "Sample_56")

growth_rates_from_FBA("Sample_68_community.pickle", "Sample_68")
add_taxonomies_to_GRs("Sample_68_fba_output_summary.csv", "Sample_68")

growth_rates_from_FBA("Sample_37_community.pickle", "Sample_37")
add_taxonomies_to_GRs("Sample_37_fba_output_summary.csv", "Sample_37")

growth_rates_from_FBA("Sample_09_community.pickle", "Sample_09")
add_taxonomies_to_GRs("Sample_09_fba_output_summary.csv", "Sample_09")

growth_rates_from_FBA("Sample_52_community.pickle", "Sample_52")
add_taxonomies_to_GRs("Sample_52_fba_output_summary.csv", "Sample_52")

growth_rates_from_FBA("Sample_95_community.pickle", "Sample_95")
add_taxonomies_to_GRs("Sample_95_fba_output_summary.csv", "Sample_95")

growth_rates_from_FBA("Sample_73_community.pickle", "Sample_73")
add_taxonomies_to_GRs("Sample_73_fba_output_summary.csv", "Sample_73")

growth_rates_from_FBA("Sample_28_community.pickle", "Sample_28")
add_taxonomies_to_GRs("Sample_28_fba_output_summary.csv", "Sample_28")

growth_rates_from_FBA("Sample_16_community.pickle", "Sample_16")
add_taxonomies_to_GRs("Sample_16_fba_output_summary.csv", "Sample_16")

growth_rates_from_FBA("Sample_67_community.pickle", "Sample_67")
add_taxonomies_to_GRs("Sample_67_fba_output_summary.csv", "Sample_67")

growth_rates_from_FBA("Sample_59_community.pickle", "Sample_59")
add_taxonomies_to_GRs("Sample_59_fba_output_summary.csv", "Sample_59")

growth_rates_from_FBA("Sample_23_community.pickle", "Sample_23")
add_taxonomies_to_GRs("Sample_23_fba_output_summary.csv", "Sample_23")

growth_rates_from_FBA("Sample_78_community.pickle", "Sample_78")
add_taxonomies_to_GRs("Sample_78_fba_output_summary.csv", "Sample_78")

growth_rates_from_FBA("Sample_46_community.pickle", "Sample_46")
add_taxonomies_to_GRs("Sample_46_fba_output_summary.csv", "Sample_46")

growth_rates_from_FBA("Sample_92_community.pickle", "Sample_92")
add_taxonomies_to_GRs("Sample_92_fba_output_summary.csv", "Sample_92")

growth_rates_from_FBA("Sample_74_community.pickle", "Sample_74")
add_taxonomies_to_GRs("Sample_74_fba_output_summary.csv", "Sample_74")

growth_rates_from_FBA("Sample_11_community.pickle", "Sample_11")
add_taxonomies_to_GRs("Sample_11_fba_output_summary.csv", "Sample_11")

growth_rates_from_FBA("Sample_30_community.pickle", "Sample_30")
add_taxonomies_to_GRs("Sample_30_fba_output_summary.csv", "Sample_30")

growth_rates_from_FBA("Sample_24_community.pickle", "Sample_24")
add_taxonomies_to_GRs("Sample_24_fba_output_summary.csv", "Sample_24")

growth_rates_from_FBA("Sample_41_community.pickle", "Sample_41")
add_taxonomies_to_GRs("Sample_41_fba_output_summary.csv", "Sample_41")

growth_rates_from_FBA("Sample_86_community.pickle", "Sample_86")
add_taxonomies_to_GRs("Sample_86_fba_output_summary.csv", "Sample_86")

growth_rates_from_FBA("Sample_05_community.pickle", "Sample_05")
add_taxonomies_to_GRs("Sample_05_fba_output_summary.csv", "Sample_05")

growth_rates_from_FBA("Sample_64_community.pickle", "Sample_64")
add_taxonomies_to_GRs("Sample_64_fba_output_summary.csv", "Sample_64")

growth_rates_from_FBA("Sample_82_community.pickle", "Sample_82")
add_taxonomies_to_GRs("Sample_82_fba_output_summary.csv", "Sample_82")

growth_rates_from_FBA("Sample_20_community.pickle", "Sample_20")
add_taxonomies_to_GRs("Sample_20_fba_output_summary.csv", "Sample_20")

growth_rates_from_FBA("Sample_45_community.pickle", "Sample_45")
add_taxonomies_to_GRs("Sample_45_fba_output_summary.csv", "Sample_45")

growth_rates_from_FBA("Sample_34_community.pickle", "Sample_34")
add_taxonomies_to_GRs("Sample_34_fba_output_summary.csv", "Sample_34")

growth_rates_from_FBA("Sample_89_community.pickle", "Sample_89")
add_taxonomies_to_GRs("Sample_89_fba_output_summary.csv", "Sample_89")

growth_rates_from_FBA("Sample_96_community.pickle", "Sample_96")
add_taxonomies_to_GRs("Sample_96_fba_output_summary.csv", "Sample_96")

growth_rates_from_FBA("Sample_70_community.pickle", "Sample_70")
add_taxonomies_to_GRs("Sample_70_fba_output_summary.csv", "Sample_70")

growth_rates_from_FBA("Sample_15_community.pickle", "Sample_15")
add_taxonomies_to_GRs("Sample_15_fba_output_summary.csv", "Sample_15")

growth_rates_from_FBA("Sample_54_community.pickle", "Sample_54")
add_taxonomies_to_GRs("Sample_54_fba_output_summary.csv", "Sample_54")

growth_rates_from_FBA("Sample_31_community.pickle", "Sample_31")
add_taxonomies_to_GRs("Sample_31_fba_output_summary.csv", "Sample_31")

growth_rates_from_FBA("Sample_10_community.pickle", "Sample_10")
add_taxonomies_to_GRs("Sample_10_fba_output_summary.csv", "Sample_10")

growth_rates_from_FBA("Sample_75_community.pickle", "Sample_75")
add_taxonomies_to_GRs("Sample_75_fba_output_summary.csv", "Sample_75")

growth_rates_from_FBA("Sample_93_community.pickle", "Sample_93")
add_taxonomies_to_GRs("Sample_93_fba_output_summary.csv", "Sample_93")

growth_rates_from_FBA("Sample_87_community.pickle", "Sample_87")
add_taxonomies_to_GRs("Sample_87_fba_output_summary.csv", "Sample_87")#

growth_rates_from_FBA("Sample_61_community.pickle", "Sample_61")
add_taxonomies_to_GRs("Sample_61_fba_output_summary.csv", "Sample_61")#

growth_rates_from_FBA("Sample_40_community.pickle", "Sample_40")
add_taxonomies_to_GRs("Sample_40_fba_output_summary.csv", "Sample_40")#

growth_rates_from_FBA("Sample_25_community.pickle", "Sample_25")
add_taxonomies_to_GRs("Sample_25_fba_output_summary.csv", "Sample_25")#

growth_rates_from_FBA("Sample_44_community.pickle", "Sample_44")
add_taxonomies_to_GRs("Sample_44_fba_output_summary.csv", "Sample_44")#

growth_rates_from_FBA("Sample_21_community.pickle", "Sample_21")
add_taxonomies_to_GRs("Sample_21_fba_output_summary.csv", "Sample_21")

growth_rates_from_FBA("Sample_83_community.pickle", "Sample_83")
add_taxonomies_to_GRs("Sample_83_fba_output_summary.csv", "Sample_83")

growth_rates_from_FBA("Sample_65_community.pickle", "Sample_65")
add_taxonomies_to_GRs("Sample_65_fba_output_summary.csv", "Sample_65")

growth_rates_from_FBA("Sample_14_community.pickle", "Sample_14")
add_taxonomies_to_GRs("Sample_14_fba_output_summary.csv", "Sample_14")

growth_rates_from_FBA("Sample_71_community.pickle", "Sample_71")
add_taxonomies_to_GRs("Sample_71_fba_output_summary.csv", "Sample_71")

growth_rates_from_FBA("Sample_50_community.pickle", "Sample_50")
add_taxonomies_to_GRs("Sample_50_fba_output_summary.csv", "Sample_50")

growth_rates_from_FBA("Sample_88_community.pickle", "Sample_88")
add_taxonomies_to_GRs("Sample_88_fba_output_summary.csv", "Sample_88")

growth_rates_from_FBA("Sample_35_community.pickle", "Sample_35")
add_taxonomies_to_GRs("Sample_35_fba_output_summary.csv", "Sample_35")





                      