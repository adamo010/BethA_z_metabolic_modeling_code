#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 09:49:41 2020

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
from micom import load_pickle
from micom.media import minimal_medium

#step 0: make a list of all the pickle files
pickle_file_list = []

for file in os.listdir():
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

#intermediate step 1: test on a single pickle file/community
sample_88_comm = load_pickle("Sample_88_community.pickle") #import the community- rate limiting step
os.system('say "Pickle import completed"')
sample_88_sol = sample_88_comm.cooperative_tradeoff(fraction = 0.5) # runs the metabolic model as tradeoffs
#okay woof,
#coop_tradeoff_sol_list.append(sample_sol) #save solution table into a list
#coop_tradeoff_sol_samples.append(key)# associate sample name with solution table.
#now, we'll work on import/export fluxes
rates = sample_88_sol.members.growth_rate.drop("medium")
medium = minimal_medium(sample_88_comm, sample_88_sol.growth_rate, min_growth = rates, exports= True) #don't know if this will work
#the above code did not, in fact, work- minimization of medium unsuccessful. Turns out we can't use the results from 
#cooperative tradeoff as additional constraints in the minimal media calculation
medium2= minimal_medium(sample_88_comm, 0.5, exports=True) #taking a long time. But, did work. 
#coop_tradeoff_med_list.append(medium)
#coop_tradeoff_sol_list.append(key)
#del sample_comm
#del sample_sol
#del rates
#del medium
#os.system('say "Another one bites the dust"')

#Intermediate step 2: write functions to get the output of FBA into files
def model_modifications(fba_output_value, fba_output_key):
    output = open("_0.5_prefba_output_summary.csv", 'w')
    results = fba_output_value.members #this is a pandas dataframe
    pd.set_option('display.max_rows', None)
    output.write(str(results))
    output.close()
    del results
    del output
    with open("_0.5_fba_output_summary.csv", "w") as csvFile:         #create a new Diet_flux_summary file with columns aligned properly etc
        writer = csv.writer(csvFile)
        with open("_0.5_prefba_output_summary.csv") as input:
            fluxreader = csv.reader(input)
            fluxreader2 = filter(None, fluxreader) #get rid of empty rows
            for row in fluxreader2:
                row = row[0].split()
                firstel = row[0]
                if len(row) < 5:
                    row = ['compartments', 'abundance', 'growth_rate', 'reactions', 'metabolites']
                writer.writerow(row)                    #create a new Diet_flux_summary file with columns aligned properly etc
    df = pd.read_csv('_0.5_fba_output_summary.csv', skiprows=1)
    df.to_csv('_0.5_fba_output_summary.csv', index=False)
    os.remove("_0.5_prefba_output_summary.csv")
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_0.5_fba_output_summary.csv"):
            dst = str(fba_output_key)+file
            os.rename(src,dst)
    del src
    del dst
    del file
    del row
    del fluxreader
    del fluxreader2
    del firstel
    return

model_modifications(sample_88_sol, "Sample_88") #tight. Don't worry about taxonomy for now.

#get the growth medium into files
def inputs_outputs(fba_medium_value, fba_output_key):
    output = open("_0.5_prefba_medium_summary.csv", 'w')
    results = fba_medium_value #this is a pandas dataframe
    pd.set_option('display.max_rows', None)
    output.write(str(results))
    output.close()
    del results
    del output
    with open("_0.5_fba_medium_summary.csv", "w") as csvFile:         #create a new Diet_flux_summary file with columns aligned properly etc
        writer = csv.writer(csvFile)
        with open("_0.5_prefba_medium_summary.csv") as input:
            fluxreader = csv.reader(input)
            fluxreader2 = filter(None, fluxreader) #get rid of empty rows
            for row in fluxreader2:
                row = row[0].split()
                firstel = row[0]
                writer.writerow(row)                    #create a new Diet_flux_summary file with columns aligned properly etc
    df = pd.read_csv('_0.5_fba_medium_summary.csv')
    df.columns=["metabolite", "flux_value"]
    df.to_csv('_0.5_fba_medium_summary.csv', index=False)
    os.remove("_0.5_prefba_medium_summary.csv")
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_0.5_fba_medium_summary.csv"):
            dst = str(fba_output_key)+file
            os.rename(src,dst)
    del src
    del dst
    del file
    del row
    del fluxreader
    del fluxreader2
    del firstel
    return

inputs_outputs(medium2, "Sample_88")

#########################################
#Step 1: select 3 pickle files to run

test_pickle_dict = {"Sample_07": "Sample_07_community.pickle", "Sample_39": "Sample_39_community.pickle", "Sample_84": "Sample_84_community.pickle"}

#this step sets up the lists we're going to populate
coop_tradeoff_sol_list = []
coop_tradeoff_sol_samples = []
coop_tradeoff_med_list = []
coop_tradeoff_med_samples = []

#Step 1: running micom
for key, value in test_pickle_dict.items():
    sample_comm = load_pickle(value) #import the community- rate limiting step
    os.system('say "Pickle import completed"')
    sample_sol = sample_comm.cooperative_tradeoff(fraction = 0.5) # runs the metabolic model as tradeoff
    coop_tradeoff_sol_list.append(sample_sol) #save solution table into a list
    coop_tradeoff_sol_samples.append(key)# associate sample name with solution table.
    #now, we'll work on import/export fluxes
    rates = sample_sol.members.growth_rate.drop("medium")
    medium = minimal_medium(sample_comm, 0.5, exports=True)
    os.system('say "Medium values calculated"')
    coop_tradeoff_med_list.append(medium)
    coop_tradeoff_med_samples.append(key)
    del sample_comm
    del sample_sol
    del rates
    del medium
os.system('say "This house is clear"')

#This step saves the MICOM coop tradeoff solution into a dictionary with its sample name as key
coop_tradeoff_dict = dict(zip(coop_tradeoff_sol_samples, coop_tradeoff_sol_list)) 
coop_tradeoff_med_dict = dict(zip(coop_tradeoff_med_samples, coop_tradeoff_med_list)) 

#this step saves the 'printed output' of the micom solution- i.e. the community growth rate (this might get funky with the table, but that's okay)
with open('10.08.20_3sample_coop_tradeoff_MICOM.csv', 'w') as csv_file:  
    writer = csv.writer(csv_file)
    for key, value in coop_tradeoff_dict.items():
       writer.writerow([key, value])

#this step saves individual species growth rates
for key, value in coop_tradeoff_dict.items():
    model_modifications(value, key)       
       
#this step saves individual species media
for key, value in coop_tradeoff_med_dict.items():
    inputs_outputs(value, key)       
       
#######################################
#phew, I don't think I can do this even for small chunks. Time to run each sample individually.
    
sample_solutions_list = []
sample_names_list = []    

def metabolic_modeling(sample_file, sample_name):
    sample_comm = load_pickle(sample_file) #import the community- rate limiting step
    sample_sol = sample_comm.cooperative_tradeoff(fraction = 0.5) # runs the metabolic model as tradeoff
    sample_solutions_list.append(sample_sol)
    sample_names_list.append(sample_name)
    rates = sample_sol.members.growth_rate.drop("medium")
    medium = minimal_medium(sample_comm, 0.5, exports=True)
    model_modifications(sample_sol, sample_name)
    inputs_outputs(medium, sample_name)
    del sample_comm
    del sample_sol
    del rates
    del medium
    os.system('say "This house is clear"')

metabolic_modeling("Sample_07_community.pickle", "Sample_07") #   
metabolic_modeling("Sample_39_community.pickle", "Sample_39") #   
metabolic_modeling("Sample_84_community.pickle", "Sample_84")
metabolic_modeling("Sample_62_community.pickle", "Sample_62")
metabolic_modeling("Sample_43_community.pickle", "Sample_43")
metabolic_modeling("Sample_18_community.pickle", "Sample_18")
metabolic_modeling("Sample_26_community.pickle", "Sample_26")
metabolic_modeling("Sample_57_community.pickle", "Sample_57")
metabolic_modeling("Sample_69_community.pickle", "Sample_69")
metabolic_modeling("Sample_32_community.pickle", "Sample_32")
metabolic_modeling("Sample_13_community.pickle", "Sample_13")
metabolic_modeling("Sample_48_community.pickle", "Sample_48")
metabolic_modeling("Sample_76_community.pickle", "Sample_76")
metabolic_modeling("Sample_90_community.pickle", "Sample_90")
metabolic_modeling("Sample_29_community.pickle", "Sample_29")
metabolic_modeling("Sample_17_community.pickle", "Sample_17")
metabolic_modeling("Sample_72_community.pickle", "Sample_72")
metabolic_modeling("Sample_94_community.pickle", "Sample_94")
metabolic_modeling("Sample_53_community.pickle", "Sample_53")
metabolic_modeling("Sample_36_community.pickle", "Sample_36")
metabolic_modeling("Sample_08_community.pickle", "Sample_08")
metabolic_modeling("Sample_79_community.pickle", "Sample_79")
metabolic_modeling("Sample_47_community.pickle", "Sample_47")
metabolic_modeling("Sample_22_community.pickle", "Sample_22")
metabolic_modeling("Sample_80_community.pickle", "Sample_80")
metabolic_modeling("Sample_66_community.pickle", "Sample_66")
metabolic_modeling("Sample_58_community.pickle", "Sample_58")
metabolic_modeling("Sample_19_community.pickle", "Sample_19")
metabolic_modeling("Sample_27_community.pickle", "Sample_27")
metabolic_modeling("Sample_42_community.pickle", "Sample_42")
metabolic_modeling("Sample_63_community.pickle", "Sample_63")
metabolic_modeling("Sample_85_community.pickle", "Sample_85")
metabolic_modeling("Sample_06_community.pickle", "Sample_06")
metabolic_modeling("Sample_38_community.pickle", "Sample_38")
metabolic_modeling("Sample_91_community.pickle", "Sample_91")
metabolic_modeling("Sample_49_community.pickle", "Sample_49")
metabolic_modeling("Sample_77_community.pickle", "Sample_77")
metabolic_modeling("Sample_12_community.pickle", "Sample_12")
metabolic_modeling("Sample_33_community.pickle", "Sample_33")
metabolic_modeling("Sample_56_community.pickle", "Sample_56")
metabolic_modeling("Sample_68_community.pickle", "Sample_68")
metabolic_modeling("Sample_37_community.pickle", "Sample_37")
metabolic_modeling("Sample_09_community.pickle", "Sample_09")
metabolic_modeling("Sample_52_community.pickle", "Sample_52")
metabolic_modeling("Sample_95_community.pickle", "Sample_95")
metabolic_modeling("Sample_73_community.pickle", "Sample_73")
metabolic_modeling("Sample_28_community.pickle", "Sample_28")
metabolic_modeling("Sample_16_community.pickle", "Sample_16")
metabolic_modeling("Sample_67_community.pickle", "Sample_67")
metabolic_modeling("Sample_59_community.pickle", "Sample_59")
metabolic_modeling("Sample_23_community.pickle", "Sample_23")
metabolic_modeling("Sample_78_community.pickle", "Sample_78")
metabolic_modeling("Sample_46_community.pickle", "Sample_46")
metabolic_modeling("Sample_92_community.pickle", "Sample_92")
metabolic_modeling("Sample_74_community.pickle", "Sample_74")
metabolic_modeling("Sample_11_community.pickle", "Sample_11")
metabolic_modeling("Sample_30_community.pickle", "Sample_30")
metabolic_modeling("Sample_24_community.pickle", "Sample_24")
metabolic_modeling("Sample_41_community.pickle", "Sample_41")
metabolic_modeling("Sample_86_community.pickle", "Sample_86")
metabolic_modeling("Sample_05_community.pickle", "Sample_05")
metabolic_modeling("Sample_64_community.pickle", "Sample_64")
metabolic_modeling("Sample_82_community.pickle", "Sample_82")
metabolic_modeling("Sample_20_community.pickle", "Sample_20")
metabolic_modeling("Sample_45_community.pickle", "Sample_45")
metabolic_modeling("Sample_34_community.pickle", "Sample_34")
metabolic_modeling("Sample_89_community.pickle", "Sample_89")
metabolic_modeling("Sample_96_community.pickle", "Sample_96")
metabolic_modeling("Sample_70_community.pickle", "Sample_70")
metabolic_modeling("Sample_15_community.pickle", "Sample_15")
metabolic_modeling("Sample_54_community.pickle", "Sample_54")
metabolic_modeling("Sample_31_community.pickle", "Sample_31")
metabolic_modeling("Sample_10_community.pickle", "Sample_10")
metabolic_modeling("Sample_75_community.pickle", "Sample_75")
metabolic_modeling("Sample_93_community.pickle", "Sample_93")
metabolic_modeling("Sample_87_community.pickle", "Sample_87")
metabolic_modeling("Sample_61_community.pickle", "Sample_61")
metabolic_modeling("Sample_40_community.pickle", "Sample_40")
metabolic_modeling("Sample_25_community.pickle", "Sample_25")
metabolic_modeling("Sample_44_community.pickle", "Sample_44")
metabolic_modeling("Sample_21_community.pickle", "Sample_21")
metabolic_modeling("Sample_83_community.pickle", "Sample_83")
metabolic_modeling("Sample_65_community.pickle", "Sample_65")
metabolic_modeling("Sample_14_community.pickle", "Sample_14")
metabolic_modeling("Sample_71_community.pickle", "Sample_71")
metabolic_modeling("Sample_50_community.pickle", "Sample_50")
metabolic_modeling("Sample_88_community.pickle", "Sample_88")
metabolic_modeling("Sample_35_community.pickle", "Sample_35")

community_growth_rate_results = dict(zip(sample_names_list, sample_solutions_list)) 
with open('10.08.20_community_GRs.csv', 'w') as csv_file:  
    writer = csv.writer(csv_file)
    for key, value in community_growth_rate_results.items():
       writer.writerow([key, value])
       
       
       
       
       