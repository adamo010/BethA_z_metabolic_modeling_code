#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 16:58:46 2022

@author: adamo010
"""
#today's goal is to see if we can clean up 02.23.22_single_sample_MICOM_to_CORDA_test_V1.py

#start with 08.26.21_V12_Hale2018_MICOM_on_MSI_fluxes_on_MMgen.py
#(I realise this isn't Hale data, it's Burns data, but the Burns MSI script is buried... somewhere)

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
import micom
from micom import load_pickle
from micom.media import minimal_medium
from itertools import chain, groupby
import cobra
from cobra.io import load_matlab_model
from cobra.medium import minimal_medium
import corda
from corda import CORDA

######################MICROBIOME/MICOM##################

##########Step 1: set wd.
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/")

##########Step 2: import file of interest
#use Patient_Blind_ID 1457; this corresponds to Sample_56/s10_S36 (tumor) and Sample_49/s09_S35 (normal)
#our pickle file: Sample56 (tumor microbiome)
imported_pickle_model = load_pickle("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Pickle_files_for_comms/Burns_pickle_files/Sample_56_community.pickle")
pickle_file_sample_name = "Sample56"

#########Step 3: calculate MM for model. 
model_only_MM =  minimal_medium(imported_pickle_model, 1.0) #old way of getting model_only_MM
#save medium information for later
medium_set3 = copy.deepcopy(model_only_MM).to_frame()
medium_set3.reset_index(inplace=True)
model_only_MM_dict = model_only_MM.to_dict() #creates a dictionary of model's coop tradeoff minimal Rxs and their fluxes

#########Step 4: import EUstandard diet
os.chdir("/Users/adamo010/Documents/microbiome_FBA/")
EUstandard = {} #editing with every iteration, so need to reset.
with open("EUstandard_fluxes.tsv", "r") as file:
    next(file)
    for line in file:                                   
        linestuff = line.strip().split()                    
        EUstandard[linestuff[0]] = float(linestuff[1])  #this creates a dictionary called 'EUstandard' containing the metabolite (key) and its amount (value) from the supplied diet file
EUstandard_edited = {k.replace("[e]", "_m"): v for k, v in EUstandard.items()}
EUstandard_edited2 = {k.replace("(e)", "_m"): v for k, v in EUstandard_edited.items()}
del file, line, linestuff, EUstandard, EUstandard_edited
#save medium information for later
EU_diet = EUstandard_edited2.items()
EU_diet_list = list(EU_diet)
medium_set1 = pd.DataFrame(EU_diet_list)
del EU_diet, EU_diet_list
#then, find shared metabolites between EU diet and metabolites that model can import
imported_pickle_model_exchange_ids = [exchange.id for exchange in imported_pickle_model.exchanges]  #list comprehension: gets all the export Rxs (exchange ids) from the  model   

#########Step 5: identify metabolites from EUstandard diet that can be imported
to_keep_from_EU_diet = {}
for key, value in EUstandard_edited2.items():
    if key in imported_pickle_model_exchange_ids:
    	to_keep_from_EU_diet[key]= value
del key, value
#save medium information for later
useable_EU_diet = to_keep_from_EU_diet.items()
useable_EU_diet_list = list(useable_EU_diet)
medium_set2 = pd.DataFrame(useable_EU_diet_list)
del useable_EU_diet, useable_EU_diet_list

#########Step 6: merge model MM and model-useable EUstd diet
keys = set(to_keep_from_EU_diet.keys()).union(model_only_MM_dict.keys())
sample_driven_diet_amounts = {k:max(to_keep_from_EU_diet.get(k,float('-inf')), model_only_MM_dict.get(k, float('-inf'))) for k in keys}
#save medium information for later
microbiome_premed = sample_driven_diet_amounts.items()
microbiome_premed_list = list(microbiome_premed)
medium_set4 = pd.DataFrame(microbiome_premed_list)
del microbiome_premed, microbiome_premed_list

#########Step 7: save all outputs
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/")
with open('_localgen_MM_from_model.csv', 'w') as f:
    for key in model_only_MM_dict.keys():
        f.write("%s,%s\n"%(key,model_only_MM_dict[key]))
del f, key
with open('_localgen_model_plus_diet_MM.csv', 'w') as f:
    for key in sample_driven_diet_amounts.keys():
    	f.write("%s,%s\n"%(key,sample_driven_diet_amounts[key]))
del f, key		

#########Step 7: save all outputs; PAUSE for now
#os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/")
#with open('_localgen_MM_from_model.csv', 'w') as f:
    #for key in model_only_MM_dict.keys():
        #f.write("%s,%s\n"%(key,model_only_MM_dict[key]))
#del f, key
#with open('_localgen_model_plus_diet_MM.csv', 'w') as f:
    #for key in sample_driven_diet_amounts.keys():
    	#f.write("%s,%s\n"%(key,sample_driven_diet_amounts[key]))
#del f, key		

#########Step 8: set minimal medium and run metabolic modeling
imported_pickle_model.medium= sample_driven_diet_amounts #used to be master_MM, but now I'm generating locally   
sample_sol = imported_pickle_model.cooperative_tradeoff(fraction=0.20, fluxes=True, pfba=False) ####################START HERE!!!!!!!
#now, extract spent medium from this.
post_coop_med = minimal_medium(imported_pickle_model, sample_sol.growth_rate, exports=True) 
#save medium information for later
medium_set5 = copy.deepcopy(post_coop_med).to_frame()
medium_set5.reset_index(inplace=True)

#########Step 9: actually calculate the change in growth medium before and after model running
#post coop med is the CHANGE in growth medium, not the growth medium itself. So we need to combine it with sample_driven_diet_amounts
#first, convert sample_driven_diet_amounts to a dataframe.
pre_micom_med = sample_driven_diet_amounts.items()
pre_micom_med_list = list(pre_micom_med)
pre_micom_med_df = pd.DataFrame(pre_micom_med_list)
#then, rename columns etc.
pre_micom_med_df.rename(columns={0: "premed_name", 1: "premed_amount"}, inplace=True)
post_coop_med_df = post_coop_med.to_frame() #convert post_coop_med (which is a series) to dataframe
post_coop_med_df.reset_index(inplace=True)
post_coop_med_df.rename(columns={"index": "postmed_name", 0: "postmed_amount"}, inplace=True)
#then, merge the dfs
spent_med_merge= pd.merge(left=pre_micom_med_df, right=post_coop_med_df, how="left", left_on="premed_name", right_on="postmed_name")
#okay, now need to double check something here. 
#By convention exports have a negative sign and imports a positive one.
#Which means exports (adding to the medium) are negative. I'll need to reverse the signs on all these, then, to get the "exported" medium
revsigns = []
for row in spent_med_merge["postmed_amount"]:
    if str(row) != "nan":
        row2 = -(row)
        revsigns.append(row2)
    else:
        revsigns.append(0)
del row, row2        
spent_med_merge["revsign_postmed_amount"] = revsigns
spent_med_merge.drop(columns={"postmed_name", "postmed_amount"}, inplace=True)
spent_med_merge["net_postmed_amount"] = spent_med_merge["premed_amount"] + spent_med_merge["revsign_postmed_amount"]
spent_med_merge.drop(columns={"revsign_postmed_amount"}, inplace=True)
#okay, now we have a df with three columns: premed_name (metabolite name), premed_amount (no micom run), and 
#net_postmed_amount (how medium has changed based on micom growth)
#let's make these into two different growth media, and run them with the human model to see if there's any difference. 
#first, let's compare EUstandard_medium with our outputs
medium_set6 = copy.deepcopy(spent_med_merge)

############Step 10: incorporate other growth media components to actually prepare spent medium
#first, find the components of model MM that are ABSENT from EUstd_diet; #this corresponds to medium_set3 vs medium_set1
medium_set7a = pd.merge(left = medium_set3, right= medium_set1, how='left', left_on="index", right_on=0, indicator=True)
#the indicator setting tells us the merge results- if the indicator returns left_only, it is only present in medium_set3 (i.e. the model_only_MM)
#these metabolites should be removed for "spent medium", as they are only present b/c they're required for models and not present in the medium
medium_set7 = medium_set7a[medium_set7a['_merge'] == "left_only"] 
medium_set7.drop(columns={"0_y", 1, "_merge"}, inplace=True)
del medium_set7a
#great. Now, we also need a list of metabolites present in the EUStandard diet that cannot be used by the model
#this corresponds to medium_set1 minus medium_set2
medium_set8a = pd.merge(left = medium_set1, right= medium_set2, how='left', left_on=0, right_on=0, indicator=True)
medium_set8 = medium_set8a[medium_set8a['_merge'] == "left_only"] 
medium_set8.drop(columns={"1_y", "_merge"}, inplace=True)
del medium_set8a
#now, we want to get spent medium, which is medium set 6-7+8
medium_set9a = pd.merge(left=medium_set6, right= medium_set7, how="left", left_on="premed_name", right_on="index", indicator=True)
#if indicator=both, that means those metabolites have been added as part of the microbial models. They should NOT be included in the "spent" medium
#for the host models. So, we're going to delete them
medium_set9b = medium_set9a[medium_set9a["_merge"] != "both"]
medium_set9b.drop(columns={"premed_amount", "index", "0_x", "_merge"}, inplace=True)
#great. Next step- add metabolites from set 8 (part of EUstd but not useable by the model. )
#can probably just stack these vertically
medium_set8.rename(columns={0: "premed_name", "1_x": "net_postmed_amount"}, inplace=True)
#now, merge for the final one!
medium_set9= pd.concat([medium_set9b, medium_set8], axis=0)
#save this medium- NOT YET
#medium_set9.to_csv("Sample_56_TD0.2_EUstd_spent_medium.csv")

######CLEAN UP!
del imported_pickle_model, pickle_file_sample_name, model_only_MM, medium_set3, model_only_MM_dict, EUstandard_edited2, sample_driven_diet_amounts,\
sample_sol, post_coop_med, medium_set5, medium_set2, pre_micom_med, pre_micom_med_list, pre_micom_med_df, post_coop_med_df, spent_med_merge, revsigns

######################HUMAN/CORDA#################

#Step 1: import human models
healthy_human_model = cobra.io.load_matlab_model("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_sample_specific_human_models/s09_S35_V3_model.mat")
tumor_human_model = cobra.io.load_matlab_model("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_sample_specific_human_models/s10_S36_V3_model.mat")

#Step 2: get a list of metabolites that the tumor model can use
max_tumor_growth= tumor_human_model.slim_optimize()
tumor_model_only_MM =  minimal_medium(tumor_human_model, max_tumor_growth) #get the model's MM
#medium_set10t = copy.deepcopy(tumor_model_only_MM).to_frame()
#medium_set10t.reset_index(inplace=True)
tumor_model_only_MM_dict = tumor_model_only_MM.to_dict() #creates a dictionary of model's minimal Rxs and their fluxes
spent_med = dict(zip(medium_set9.premed_name, medium_set9.net_postmed_amount)) #make a diet dictionary of metabolites and their amounts
tumor_model_exchange_ids = [exchange.id for exchange in tumor_human_model.exchanges]  #list comprehension: gets all the export Rxs (exchange ids) from the  model   
tumor_to_keep_from_spentmed = {} #initiate a dictionary of metabs/amounts to keep from highfiber dict
for key, value in spent_med.items():
    if key in tumor_model_exchange_ids:
        tumor_to_keep_from_spentmed[key]= value
del key, value
#save metabolites from spent med which are useable by model. NOT CURRENTLY USED
#tumor_useable_spent_med = tumor_to_keep_from_spentmed.items()
#tumor_useable_spent_med_list = list(tumor_useable_spent_med)
#medium_set11t = pd.DataFrame(tumor_useable_spent_med_list)
#del tumor_useable_spent_med, tumor_useable_spent_med_list

#Step 3: get a diet that the sample can grow on
keys = set(tumor_to_keep_from_spentmed.keys()).union(tumor_model_only_MM_dict.keys()) #identify metabolites from diet that can be imported by model
tumor_driven_diet_amounts = {k:max(tumor_to_keep_from_spentmed.get(k,float('-inf')), tumor_model_only_MM_dict.get(k, float('-inf'))) for k in keys} #merge these dicts
#save metabolites for later: NOT CURRENTLY USED
#tumor_premed = tumor_driven_diet_amounts.items()
#tumor_premed_list = list(tumor_premed)
#medium_set12t = pd.DataFrame(tumor_premed_list)
#del tumor_premed, tumor_premed_list

#Step 4: Run the model!
tumor_human_model.medium= tumor_driven_diet_amounts #ding   
human_tumor_sample_sol = tumor_human_model.optimize() #rate limiting step; might take a while.

#Step 5: save the outputs
#first, fluxes
tumor_fluxes = human_tumor_sample_sol.fluxes
tumor_fluxes.to_csv("s10_S36_tumor_cell_tumor_microbiome_EUstd_fluxes.csv")
#them, growth rate/biomass
tumor_GR = str(human_tumor_sample_sol)
with open("s10_S36_tumor_cell_tumor_microbiome_EUstd_GR.csv", "w") as text_file:
    text_file.write(tumor_GR)
del text_file
#then, inputs/outputs
tumor_post_growth_med = minimal_medium(tumor_human_model, exports=True)
tumor_post_growth_med.to_csv("s10_S36_tumor_cell_tumor_microbiome_EUstd_uptakes_and_secretions.csv")

#Step 6: get a list of metabolites that the healthy model can use
max_healthy_growth= healthy_human_model.slim_optimize()
healthy_model_only_MM =  minimal_medium(healthy_human_model, max_healthy_growth) #get the model's MM
#medium_set10h = copy.deepcopy(healthy_model_only_MM).to_frame()
#medium_set10h.reset_index(inplace=True)
healthy_model_only_MM_dict = healthy_model_only_MM.to_dict() #creates a dictionary of model's minimal Rxs and their fluxes
spent_med = dict(zip(medium_set9.premed_name, medium_set9.net_postmed_amount)) #make a diet dictionary of metabolites and their amounts
healthy_model_exchange_ids = [exchange.id for exchange in healthy_human_model.exchanges]  #list comprehension: gets all the export Rxs (exchange ids) from the  model   
healthy_to_keep_from_spentmed = {} #initiate a dictionary of metabs/amounts to keep from highfiber dict
for key, value in spent_med.items():
    if key in healthy_model_exchange_ids:
        healthy_to_keep_from_spentmed[key]= value
del key, value
#save metabolites from spent med which are useable by model. NOT CURRENTLY USED
#healthy_useable_spent_med = healthy_to_keep_from_spentmed.items()
#healthy_useable_spent_med_list = list(healthy_useable_spent_med)
#medium_set11h = pd.DataFrame(healthy_useable_spent_med_list)
#del healthy_useable_spent_med, healthy_useable_spent_med_list

#Step 7: get a diet that the sample can grow on
keys = set(healthy_to_keep_from_spentmed.keys()).union(healthy_model_only_MM_dict.keys()) #identify metabolites from diet that can be imported by model
healthy_driven_diet_amounts = {k:max(healthy_to_keep_from_spentmed.get(k,float('-inf')), healthy_model_only_MM_dict.get(k, float('-inf'))) for k in keys} #merge these dicts
#save metabolites:NOT CURRENTLY USED
#healthy_premed = healthy_driven_diet_amounts.items()
#healthy_premed_list = list(healthy_premed)
#medium_set12h = pd.DataFrame(healthy_premed_list)
#del healthy_premed, healthy_premed_list

#Step 8: Run the model!
healthy_human_model.medium= healthy_driven_diet_amounts #ding   
human_healthy_sample_sol = healthy_human_model.optimize() #rate limiting step; might take a while.

#Step 9: save the outputs
#first, fluxes
healthy_fluxes = human_healthy_sample_sol.fluxes
healthy_fluxes.to_csv("s57_S76_healthy_cell_tumor_microbiome_EUstd_fluxes.csv")
#them, growth rate/biomass
healthy_GR = str(human_healthy_sample_sol)
with open("s57_S76_healthy_cell_tumor_microbiome_EUstd_GR.csv", "w") as text_file:
    text_file.write(healthy_GR)
del text_file
#then, inputs/outputs
healthy_post_growth_med = minimal_medium(healthy_human_model, exports=True)
healthy_post_growth_med.to_csv("s57_S76_healthy_cell_tumor_microbiome_EUstd_uptakes_and_secretions.csv")
#as a reminder, exports are negative and imports are positive




