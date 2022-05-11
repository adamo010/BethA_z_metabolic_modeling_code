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
import pickle
import sys #for specifying input files externally
from micom import load_pickle
from micom.media import minimal_medium
from itertools import chain, groupby

#order of operations:
#1: write a function to import pickle files and run micom, generating temporary output files (including community and species specific growth rates)
#2: send those files to a specific subfolder 
#3: repeat.

#we start in Combining_Hale_models

def main():
    script = sys.argv[0]
    pickle_file = sys.argv[1]
    pickle_file_sample_prename = str(sys.argv[1])
    pickle_file_sample_name = pickle_file_sample_prename.replace("_community.pickle", "")
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Niccolai_models/Niccolai2020_pickle_files/")
    imported_pickle_model = load_pickle(pickle_file)
    ###Step 1: calculate MM for model. 
    model_only_MM =  minimal_medium(imported_pickle_model, 1.0) #old way of getting model_only_MM; go back to this in V16.
    model_only_MM_dict = model_only_MM.to_dict() #creates a dictionary of model's coop tradeoff minimal Rxs and their fluxes
    ###Step 2: import EUstandard diet
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Niccolai_models/")
    EUstandard = {} #editing with every iteration, so need to reset.
    with open("EUstandard_fluxes.tsv", "r") as file:
        next(file)
        for line in file:                                   
            linestuff = line.strip().split()                    
            EUstandard[linestuff[0]] = float(linestuff[1])  #this creates a dictionary called 'EUstandard' containing the metabolite (key) and its amount (value) from the supplied diet file
    EUstandard_edited = {k.replace("[e]", "_m"): v for k, v in EUstandard.items()}
    EUstandard_edited2 = {k.replace("(e)", "_m"): v for k, v in EUstandard_edited.items()}
    del file, line, linestuff
    #then, find shared metabolites between EU diet and metabolites that model can import
    imported_pickle_model_exchange_ids = [exchange.id for exchange in imported_pickle_model.exchanges]  #list comprehension: gets all the export Rxs (exchange ids) from the  model   
    ###Step 3: identify metabolites from EUstandard diet that can be imported
    to_keep_from_EU_diet = {}
    for key, value in EUstandard_edited2.items():
    	if key in imported_pickle_model_exchange_ids:
    		to_keep_from_EU_diet[key]= value
    del key, value
    ###Step 4 NEW: straight up merge these dicts
    keys = set(to_keep_from_EU_diet.keys()).union(model_only_MM_dict.keys())
    sample_driven_diet_amounts = {k:max(to_keep_from_EU_diet.get(k,float('-inf')), model_only_MM_dict.get(k, float('-inf'))) for k in keys}
    #step 6: save all outputs
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Niccolai_models/Niccolai_MICOM_output/")
    with open('_V01_Niccolai_MSIgen_MM_from_model.csv', 'w') as f:
        for key in model_only_MM_dict.keys():
            f.write("%s,%s\n"%(key,model_only_MM_dict[key]))
    del f, key
    with open('_V01_Niccolai_MSIgen_model_plus_diet_MM.csv', 'w') as f:
    	for key in sample_driven_diet_amounts.keys():
    		f.write("%s,%s\n"%(key,sample_driven_diet_amounts[key]))
    del f, key		
    ###Step 7: set minimal medium and run metabolic modeling
    imported_pickle_model.medium= sample_driven_diet_amounts #used to be master_MM, but now I'm generating locally   
    sample_sol = imported_pickle_model.cooperative_tradeoff(fraction=0.25, fluxes=True, pfba=False) 
    rates = sample_sol.members.growth_rate.drop("medium")
    fluxes = sample_sol.fluxes
    post_coop_med = minimal_medium(imported_pickle_model, sample_sol.growth_rate, exports=True) 
    comm_GR = str(sample_sol)
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Niccolai_models/Niccolai_MICOM_output/")
    rates.to_csv("_V01_Niccolai_indiv_spp_GRs.csv")
    fluxes.to_csv("_V01_Niccolai_fluxes.csv")
    if post_coop_med is None:
    	output = "This sample has no post coop med"
    	with open("_V01_Niccolai_input_output.csv", "w") as text_file:
    		text_file.write(output)    
    else:
    	post_coop_med.to_csv("_V01_Niccolai_input_output.csv")
    with open("_V01_Niccolai_comm_GRs.csv", "w") as text_file:
        text_file.write(comm_GR)
    #now, rename all files
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_V01_Niccolai_MSIgen_MM_from_model.csv"):
            dst=str(pickle_file_sample_name)+file
            os.rename(src,dst)
        elif fnmatch.fnmatch(file, "_V01_Niccolai_MSIgen_model_plus_diet_MM.csv"):
        	dst=str(pickle_file_sample_name)+file
        	os.rename(src,dst)
        elif fnmatch.fnmatch(file, "_V01_Niccolai_indiv_spp_GRs.csv"):
            dst=str(pickle_file_sample_name)+file
            os.rename(src,dst) 
        elif fnmatch.fnmatch(file, "_V01_Niccolai_fluxes.csv"):
            dst=str(pickle_file_sample_name)+file
            os.rename(src,dst)    
        elif fnmatch.fnmatch(file, "_V01_Niccolai_input_output.csv"):
        	dst=str(pickle_file_sample_name)+file
        	os.rename(src,dst)
        elif fnmatch.fnmatch(file, "_V01_Niccolai_comm_GRs.csv"):
        	dst=str(pickle_file_sample_name)+file
        	os.rename(src,dst)	   	            	       			
    return

main()
