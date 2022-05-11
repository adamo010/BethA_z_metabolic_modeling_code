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
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/Hale2018_pickle_files/")
    imported_pickle_model = load_pickle(pickle_file)
    ###Step 1: calculate MM for model. 
    model_only_MM =  minimal_medium(imported_pickle_model, 1.0) #old way of getting model_only_MM; go back to this in V16.
    model_only_MM_dict = model_only_MM.to_dict() #creates a dictionary of model's coop tradeoff minimal Rxs and their fluxes
    ###Step 2: import EUstandard diet
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/")
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
    ###Step 4: create a dictionary of the metabolites and their amounts that ARE in the model_MM list, but are NOT in the diet input
    #missing_food = {}
    #for key, value in model_only_MM_dict.items():
        #if key in to_keep_from_EU_diet:
            #pass
        #else:
            #missing_food[key] = value
    #del key, value
    ###Step 5: merge to_keep_from_EU_diet (present in EUstd and importable into model) with missing_food (found in model MM but not in EUstd)
    # = dict(**to_keep_from_EU_diet, **missing_food)
    ###Step 6: match amounts from diet_input_rounded to amounts in a minimal medium for 1g biomass.
    #x = model_only_MM_dict
    #y = sample_driven_diet_amounts
    #diet_input_min_growth = {key: max(x[key], value) if key in x.keys() else y[key] for key, value in y.items()}
    #del x, y
    #Step 5: intermediate step: figure out if models will run AT ALL????
    #premed_sample_sol = imported_pickle_model.cooperative_tradeoff(fraction=0.20, fluxes=True, pfba=True)
    #premed_rates = premed_sample_sol.members.growth_rate.drop("medium")
    #premed_fluxes= premed_sample_sol.fluxes
    #premed_post_coop_med = minimal_medium(imported_pickle_model, premed_sample_sol.growth_rate, min_growth=premed_rates, exports=True)
    #premed_comm_GR = str(premed_sample_sol)
    #step 6: save all outputs
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/Hale_MICOM_output/")
    with open('_V9_Hale_MSIgen_MM_from_model.csv', 'w') as f:
        for key in model_only_MM_dict.keys():
            f.write("%s,%s\n"%(key,model_only_MM_dict[key]))
    del f, key
    with open('_V9_Hale_MSIgen_model_plus_diet_MM.csv', 'w') as f:
    	for key in sample_driven_diet_amounts.keys():
    		f.write("%s,%s\n"%(key,sample_driven_diet_amounts[key]))
    del f, key		
    #premed_rates.to_csv("_V9_Hale_indiv_spp_GRs_noTD.csv")
    #premed_fluxes.to_csv("_V9_Hale_fluxes_noTD.csv")
    #if premed_post_coop_med is None:
    	#output = "This sample has no post coop med"
    	#with open("_V9_Hale_input_output_noTD.csv", "w") as text_file:
    		#text_file.write(output)    
    #else:
    	#premed_post_coop_med.to_csv("_V9_Hale_input_output_noTD.csv")
    #with open("_V9_Hale_comm_GRs_noTD.csv", "w") as text_file:
        #text_file.write(premed_comm_GR)     
    ###Step 7: set minimal medium and run metabolic modeling
    imported_pickle_model.medium= sample_driven_diet_amounts #used to be master_MM, but now I'm generating locally   
    sample_sol = imported_pickle_model.cooperative_tradeoff(fraction=0.20, fluxes=True, pfba=False) 
    rates = sample_sol.members.growth_rate.drop("medium")
    fluxes = sample_sol.fluxes
    post_coop_med = minimal_medium(imported_pickle_model, sample_sol.growth_rate, min_growth=rates, exports=True) 
    comm_GR = str(sample_sol)
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/Hale_MICOM_output/")
    rates.to_csv("_V9_Hale_indiv_spp_GRs.csv")
    fluxes.to_csv("_V9_Hale_fluxes.csv")
    if post_coop_med is None:
    	output = "This sample has no post coop med"
    	with open("_V9_Hale_input_output.csv", "w") as text_file:
    		text_file.write(output)    
    else:
    	post_coop_med.to_csv("_V9_Hale_input_output.csv")
    with open("_V9_Hale_comm_GRs.csv", "w") as text_file:
        text_file.write(comm_GR)
    #now, rename all files
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_V9_Hale_MSIgen_MM_from_model.csv"):
            dst=str(pickle_file_sample_name)+file
            os.rename(src,dst)
        elif fnmatch.fnmatch(file, "_V9_Hale_MSIgen_model_plus_diet_MM.csv"):
        	dst=str(pickle_file_sample_name)+file
        	os.rename(src,dst)
        #elif fnmatch.fnmatch(file, "_V9_Hale_indiv_spp_GRs_noTD.csv"):
            #dst=str(pickle_file_sample_name)+file
            #os.rename(src,dst) 
        #elif fnmatch.fnmatch(file, "_V9_Hale_fluxes_noTD.csv"):
            #dst=str(pickle_file_sample_name)+file
            #os.rename(src,dst)
        #elif fnmatch.fnmatch(file, "_V9_Hale_input_output_noTD.csv"):
        	#dst=str(pickle_file_sample_name)+file
        	#os.rename(src,dst)   	            
        #elif fnmatch.fnmatch(file, "_V9_Hale_comm_GRs_noTD.csv"):
        	#dst=str(pickle_file_sample_name)+file
        	#os.rename(src,dst)       			
        elif fnmatch.fnmatch(file, "_V9_Hale_indiv_spp_GRs.csv"):
            dst=str(pickle_file_sample_name)+file
            os.rename(src,dst) 
        elif fnmatch.fnmatch(file, "_V9_Hale_fluxes.csv"):
            dst=str(pickle_file_sample_name)+file
            os.rename(src,dst)    
        elif fnmatch.fnmatch(file, "_V9_Hale_input_output.csv"):
        	dst=str(pickle_file_sample_name)+file
        	os.rename(src,dst)
        elif fnmatch.fnmatch(file, "_V9_Hale_comm_GRs.csv"):
        	dst=str(pickle_file_sample_name)+file
        	os.rename(src,dst)	   	            	       			
    return

main()
