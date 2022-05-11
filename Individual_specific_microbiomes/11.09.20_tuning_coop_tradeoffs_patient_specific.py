#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 20:45:05 2020

@author: adamo010
"""

######I HAVE RETURNED

#so, one of the things I'd like to try now is running through a few different tradeoff parameters and seeing what impact this has on
#growth rates, who grows when, etc.

#but, if I have learned one thing, it's that this stupid process is insanely long.
#so I will work with four patients only (8 samples)

#see 10.08.20_local_patient_specific_community_GRs_clean.py for borrowing this code.

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

#so, to get a range of tradeoff parameters from 0.1 to 1.0, use the following code:
#sols = com.cooperative_tradeoff(fraction=np.arange(0.1, 1.01, 0.1))

#see MICOM_healthy_vs_CRC to figure out how to put these data into a useful dataframe

#will probably have to go sample by sample, and run MICOM on a list of samples
#ugh, each sample is going to have 10 different output files. How in the hell am I going to iteratively name them? A FUNCTION, no doubt. Gross.

#steps:
#1) identify some good patients and their respective tumor pairs
#2) Edit the metabolic modeling function in 10.08.20_local_patient_specific_community_GRs_clean.py to account for multiple tradeoff parameters
#3) Edit the other functions to account for multiple tradeoff parameters and spit out multiple metabolite and species-specific GR output files
#4) ADD IN FLUXES because I want to burn out my computer.

#the files we want to use are as follows:

sample_sub_dict = {"Sample_88": "Sample_88_community.pickle", "Sample_27": "Sample_27_community.pickle",
                   "Sample_43": "Sample_43_community.pickle", "Sample_26": "Sample_26_community.pickle", 
                   "Sample_34": "Sample_34_community.pickle", "Sample_64": "Sample_64_community.pickle",
                   "Sample_36": "Sample_36_community.pickle", "Sample_50": "Sample_50_community.pickle"}

#you know what, it might be easier to just run this one sample at a time until I get a sense of what the output is going to look like

sample_solutions_list = []
sample_names_list = []    

%time Sample_88_comm = load_pickle("Pickle_files_for_comms/Sample_88_community.pickle")
#psa this took 6min 17s
%time sample_sol = Sample_88_comm.cooperative_tradeoff(fraction=np.arange(0.1, 1.01, 0.1))
sample_solutions_list.append(str(sample_sol))
sample_names_list.append("Sample_88")

rates = sample_sol.solution.apply(lambda x: x.members.growth_rate) #fast
#aha, well, this is something, at least. Rates is a pandas dataframe where each row is a compartment(correxponding to a
#tradeoff parameter) and each column is an OTU id. The values in the table are individual taxa growth rates.
#want to save this and my overall community GRs. 

#couple of things to do with the dataframe first. Edit the first column header to read "tradeoff_value", 
#and edit the column values to be 0.1-1.0 instead of 0-9
#also, running into that issue where most of the columns are "..." instead of data.

#change row names (note that index=0 means tradeoff=1.0, index=9 means tradeoff=0.1)
rates.index = ["1.0", "0.9", "0.8", "0.7", "0.6", "0.5", "0.4", "0.3", "0.2", "0.1"]
rates=rates.drop(["medium"], axis=1)
#this is way faster than the other thing I was doing. 

rates.to_csv("Sample_88_fba_output_summary_cooptradeoff.csv")

#for file in os.listdir():
  #  src=file
   # if fnmatch.fnmatch(file, "_0.5_fba_output_summary.csv"):
      #  dst = str(fba_output_key)+file
       # os.rename(src,dst)            
#add this later, if we need to put this crap into a function.
comm_GRs = []

for row in sample_sol["solution"]:
    rowstring = str(row)
    rowstring= rowstring.split()
    print(rowstring)
    comm_GR = rowstring[1]
    comm_GRs.append(comm_GR)

community_GRs_df = copy.deepcopy(sample_sol)
community_GRs_df["solution"] = comm_GRs    
community_GRs_df.to_csv("Sample_88_fba_community_GRs_cooptradeoff.csv")
    
#######################Excellent!
#now let's run this again with the paired sample and BUCKLE THE FUCKLE UPPLE for another 24 hour run.

%time Sample_27_comm = load_pickle("Pickle_files_for_comms/Sample_27_community.pickle")
#psa this took 5min 19s
%time sample_sol = Sample_27_comm.cooperative_tradeoff(fraction=np.arange(0.1, 1.01, 0.1))
sample_solutions_list.append(str(sample_sol))
sample_names_list.append("Sample_27")

rates = sample_sol.solution.apply(lambda x: x.members.growth_rate) #fast

#change row names (note that index=0 means tradeoff=1.0, index=9 means tradeoff=0.1)
rates.index = ["1.0", "0.9", "0.8", "0.7", "0.6", "0.5", "0.4", "0.3", "0.2", "0.1"]
rates=rates.drop(["medium"], axis=1)
#this is way faster than the other thing I was doing. 

rates.to_csv("Sample_27_fba_output_summary_cooptradeoff.csv")

comm_GRs = []

for row in sample_sol["solution"]:
    rowstring = str(row)
    rowstring= rowstring.split()
    print(rowstring)
    comm_GR = rowstring[1]
    comm_GRs.append(comm_GR)

community_GRs_df = copy.deepcopy(sample_sol)
community_GRs_df["solution"] = comm_GRs    
community_GRs_df.to_csv("Sample_27_fba_community_GRs_cooptradeoff.csv")









#in retrospect, I don't actually know how we would run every cooperative tradeoff parameter for media. 
#probably some kind of foreloop. 

#okay, here is our model running function:

def metabolic_modeling(sample_file, sample_name):
    sample_comm = load_pickle(sample_file) #import the community- rate limiting step
    sample_sol = sample_comm.cooperative_tradeoff(fraction=np.arange(0.1, 1.01, 0.1)) # runs the metabolic model as tradeoff OOF
    
    sample_solutions_list.append(str(sample_sol))
    sample_names_list.append(sample_name)
    rates = sample_sol.members.growth_rate.drop("medium")
    medium = minimal_medium(sample_comm, 0.5, exports=True)
    model_modifications(sample_sol, sample_name)
    inputs_outputs(medium, sample_name)
    del sample_comm
    del sample_sol
    del rates
    del medium
    os.system('say "Moving on, we are all done here."')

