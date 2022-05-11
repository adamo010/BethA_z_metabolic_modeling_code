#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 12:40:15 2022

@author: adamo010
"""

#this is a cleaned version of 02.28.22_few_sample_MICOM_to_CORDA_test_V3.py

import gurobipy
import scipy
import numpy as np
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

#lines with delloop should be removed in MSI version; lines with addloop should be added in MSI version

###############STEP 0: file lists and names################
#have created a new prime metadata file that contains all necessary file names and metadata for this project
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/")
all_metadata= pd.read_csv("Metadata_for_human_and_microbial_communities_Burns_data.csv")
#because we're using a subset (and running the whole thing on MSI later), let's trim up the metadata a bit.
PBIDS = [1157,1178,1183] #delloop
sub_metadata = all_metadata[all_metadata["Patient_Blind_ID"].isin(PBIDS)] #delloop
del PBIDS #delloop
#will need to think about how to split this up into dictionaries for running. BUT, it is most certainly doable. 
#need subdicts for each type of model (tumor/normal)
sub_metadata_tumor = sub_metadata[sub_metadata["Description"]=="tumor"] #delloop
sub_metadata_normal = sub_metadata[sub_metadata["Description"]=="normal"] #delloop
#addloop#metadata_tumor= all_metadata[all_metadata["Description"]=="tumor"]
#addloop#metadata_healthy= all_metadata[all_metadata["Description"]=="normal"]

###############STEP 1: import and ready the EUstandard diet info#########
#note that this is out of order from previous file, but I think that's okay
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

################STEP 2: calculating MM for model. our first loop. 
#prepare dictionaries for looping
normal_microbiome_model_dict = dict(zip(sub_metadata_normal.SampleID, sub_metadata_normal.microbiome_model_file))
tumor_microbiome_model_dict = dict(zip(sub_metadata_tumor.SampleID, sub_metadata_tumor.microbiome_model_file))
sample_to_PBID_dict = dict(zip(sub_metadata.SampleID, sub_metadata.Patient_Blind_ID)) #key is sampleid (unique), value is PBID (not unique)
#move to the directory where we'll dump our output (spent med) files
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/")

################STEP 2a: normal microbiomes
for key, value in normal_microbiome_model_dict.items():
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Pickle_files_for_comms/Burns_pickle_files/") #Move to pickle file directory
    imported_pickle_model = load_pickle(value) #import our pickle file
    pickle_file_sample_name = str(key) #note the sample name
    model_only_MM =  minimal_medium(imported_pickle_model, 1.0) #get model_only_MM
    medium_set3 = copy.deepcopy(model_only_MM).to_frame() #save medium information for later
    medium_set3.reset_index(inplace=True) #clean up the dataframe
    model_only_MM_dict = model_only_MM.to_dict() #creates a dictionary of model's coop tradeoff minimal Rxs and their fluxes
    imported_pickle_model_exchange_ids = [exchange.id for exchange in imported_pickle_model.exchanges] #find shared metabolites between EU diet and metabolites that model can import
    to_keep_from_EU_diet = {} #create a new dict
    for key, value in EUstandard_edited2.items():
        if key in imported_pickle_model_exchange_ids:
            to_keep_from_EU_diet[key]= value
    del key, value #only keep items from EUstandard that match model exchange ids (i.e. can be imported by model) 
    useable_EU_diet = to_keep_from_EU_diet.items() #save medium information for later
    useable_EU_diet_list = list(useable_EU_diet)
    medium_set2 = pd.DataFrame(useable_EU_diet_list)
    del useable_EU_diet, useable_EU_diet_list
    keys = set(to_keep_from_EU_diet.keys()).union(model_only_MM_dict.keys()) #create a union of model MM and model-useable EUstd diet keys
    sample_driven_diet_amounts = {k:max(to_keep_from_EU_diet.get(k,float('-inf')), model_only_MM_dict.get(k, float('-inf'))) for k in keys} #merge dicts
    microbiome_premed = sample_driven_diet_amounts.items() #save medium information for later
    microbiome_premed_list = list(microbiome_premed)
    medium_set4 = pd.DataFrame(microbiome_premed_list)
    del microbiome_premed, microbiome_premed_list, keys
    ########Save any important medium output files
    imported_pickle_model.medium= sample_driven_diet_amounts #set sample_driven_diet_amounts to model's growth medium
    sample_sol = imported_pickle_model.cooperative_tradeoff(fraction=0.20, fluxes=True, pfba=False) #run FBA with TD=0.2
    post_coop_med = minimal_medium(imported_pickle_model, sample_sol.growth_rate, exports=True)#now, extract spent medium from this.
    medium_set5 = copy.deepcopy(post_coop_med).to_frame() #save medium information for later
    medium_set5.reset_index(inplace=True) #clean up the dataframe
    pre_micom_med = sample_driven_diet_amounts.items() #convert sample_driven_diet_amounts to a dataframe.
    pre_micom_med_list = list(pre_micom_med)
    pre_micom_med_df = pd.DataFrame(pre_micom_med_list)
    pre_micom_med_df.rename(columns={0: "premed_name", 1: "premed_amount"}, inplace=True) #then, rename columns etc.
    post_coop_med_df = post_coop_med.to_frame() #convert post_coop_med (which is a series) to dataframe
    post_coop_med_df.reset_index(inplace=True)
    post_coop_med_df.rename(columns={"index": "postmed_name", 0: "postmed_amount"}, inplace=True)
    spent_med_merge= pd.merge(left=pre_micom_med_df, right=post_coop_med_df, how="left", left_on="premed_name", right_on="postmed_name") #then, merge the dfs
    #Which means exports (adding to the medium) are negative. I'll need to reverse the signs on all these, then, to get the "exported" medium
    revsigns = [] #reverse signs to get correct values for exported medium
    for row in spent_med_merge["postmed_amount"]:
        if str(row) != "nan":
            row2 = -(row)
            revsigns.append(row2)
        else:
            revsigns.append(0)
    del row, row2       
    spent_med_merge["revsign_postmed_amount"] = revsigns #set correctly signed values as column
    spent_med_merge.drop(columns={"postmed_name", "postmed_amount"}, inplace=True)
    spent_med_merge["net_postmed_amount"] = spent_med_merge["premed_amount"] + spent_med_merge["revsign_postmed_amount"]
    spent_med_merge.drop(columns={"revsign_postmed_amount"}, inplace=True)
    medium_set6 = copy.deepcopy(spent_med_merge) #save for later
    #now, we need to incorporate other growth media components to actually prepare spent medium
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
    medium_set8.rename(columns={0: "premed_name", "1_x": "net_postmed_amount"}, inplace=True) #stack EUstd diets to merge vertically
    medium_set9= pd.concat([medium_set9b, medium_set8], axis=0) #now, merge for the final one!
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_results/")
    medium_set9.to_csv("_TD0.2_EUstd_spent_medium.csv")  #save this medium
    #pull the PBID associated with the sampleID and use it to create a unique sample name
    pbid = sample_to_PBID_dict.get(pickle_file_sample_name) #ectract PBID (value) from sample_to_PBID_dict
    unique_file_name=  "PBID_"+ str(pbid)+"_normal_microbiome"
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_TD0.2_EUstd_spent_medium.csv"):
            dst=str(unique_file_name)+file
            os.rename(src,dst)
######CLEAN UP!
del imported_pickle_model, pickle_file_sample_name, model_only_MM, medium_set3, model_only_MM_dict, sample_driven_diet_amounts,\
        sample_sol, post_coop_med, medium_set5, medium_set2, pre_micom_med, pre_micom_med_list, pre_micom_med_df, post_coop_med_df, spent_med_merge, revsigns

################STEP 2b: tumor microbiomes
for key, value in tumor_microbiome_model_dict.items():
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Pickle_files_for_comms/Burns_pickle_files/") #Move to pickle file directory
    imported_pickle_model = load_pickle(value) #import our pickle file
    pickle_file_sample_name = str(key) #note the sample name
    model_only_MM =  minimal_medium(imported_pickle_model, 1.0) #get model_only_MM
    medium_set3 = copy.deepcopy(model_only_MM).to_frame() #save medium information for later
    medium_set3.reset_index(inplace=True) #clean up the dataframe
    model_only_MM_dict = model_only_MM.to_dict() #creates a dictionary of model's coop tradeoff minimal Rxs and their fluxes
    imported_pickle_model_exchange_ids = [exchange.id for exchange in imported_pickle_model.exchanges] #find shared metabolites between EU diet and metabolites that model can import
    to_keep_from_EU_diet = {} #create a new dict
    for key, value in EUstandard_edited2.items():
        if key in imported_pickle_model_exchange_ids:
            to_keep_from_EU_diet[key]= value
    del key, value #only keep items from EUstandard that match model exchange ids (i.e. can be imported by model) 
    useable_EU_diet = to_keep_from_EU_diet.items() #save medium information for later
    useable_EU_diet_list = list(useable_EU_diet)
    medium_set2 = pd.DataFrame(useable_EU_diet_list)
    del useable_EU_diet, useable_EU_diet_list
    keys = set(to_keep_from_EU_diet.keys()).union(model_only_MM_dict.keys()) #create a union of model MM and model-useable EUstd diet keys
    sample_driven_diet_amounts = {k:max(to_keep_from_EU_diet.get(k,float('-inf')), model_only_MM_dict.get(k, float('-inf'))) for k in keys} #merge dicts
    microbiome_premed = sample_driven_diet_amounts.items() #save medium information for later
    microbiome_premed_list = list(microbiome_premed)
    medium_set4 = pd.DataFrame(microbiome_premed_list)
    del microbiome_premed, microbiome_premed_list, keys
    ########Save any important medium output files
    imported_pickle_model.medium= sample_driven_diet_amounts #set sample_driven_diet_amounts to model's growth medium
    sample_sol = imported_pickle_model.cooperative_tradeoff(fraction=0.20, fluxes=True, pfba=False) #run FBA with TD=0.2
    post_coop_med = minimal_medium(imported_pickle_model, sample_sol.growth_rate, exports=True)#now, extract spent medium from this.
    medium_set5 = copy.deepcopy(post_coop_med).to_frame() #save medium information for later
    medium_set5.reset_index(inplace=True) #clean up the dataframe
    pre_micom_med = sample_driven_diet_amounts.items() #convert sample_driven_diet_amounts to a dataframe.
    pre_micom_med_list = list(pre_micom_med)
    pre_micom_med_df = pd.DataFrame(pre_micom_med_list)
    pre_micom_med_df.rename(columns={0: "premed_name", 1: "premed_amount"}, inplace=True) #then, rename columns etc.
    post_coop_med_df = post_coop_med.to_frame() #convert post_coop_med (which is a series) to dataframe
    post_coop_med_df.reset_index(inplace=True)
    post_coop_med_df.rename(columns={"index": "postmed_name", 0: "postmed_amount"}, inplace=True)
    spent_med_merge= pd.merge(left=pre_micom_med_df, right=post_coop_med_df, how="left", left_on="premed_name", right_on="postmed_name") #then, merge the dfs
    #Which means exports (adding to the medium) are negative. I'll need to reverse the signs on all these, then, to get the "exported" medium
    revsigns = [] #reverse signs to get correct values for exported medium
    for row in spent_med_merge["postmed_amount"]:
        if str(row) != "nan":
            row2 = -(row)
            revsigns.append(row2)
        else:
            revsigns.append(0)
    del row, row2       
    spent_med_merge["revsign_postmed_amount"] = revsigns #set correctly signed values as column
    spent_med_merge.drop(columns={"postmed_name", "postmed_amount"}, inplace=True)
    spent_med_merge["net_postmed_amount"] = spent_med_merge["premed_amount"] + spent_med_merge["revsign_postmed_amount"]
    spent_med_merge.drop(columns={"revsign_postmed_amount"}, inplace=True)
    medium_set6 = copy.deepcopy(spent_med_merge) #save for later
    #now, we need to incorporate other growth media components to actually prepare spent medium
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
    medium_set8.rename(columns={0: "premed_name", "1_x": "net_postmed_amount"}, inplace=True) #stack EUstd diets to merge vertically
    medium_set9= pd.concat([medium_set9b, medium_set8], axis=0) #now, merge for the final one!
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_results/")
    medium_set9.to_csv("_TD0.2_EUstd_spent_medium.csv")  #save this medium
    #pull the PBID associated with the sampleID and use it to create a unique sample name
    pbid = sample_to_PBID_dict.get(pickle_file_sample_name) #ectract PBID (value) from sample_to_PBID_dict
    unique_file_name=  "PBID_"+ str(pbid)+"_tumor_microbiome"
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_TD0.2_EUstd_spent_medium.csv"):
            dst=str(unique_file_name)+file
            os.rename(src,dst)
######CLEAN UP!
del imported_pickle_model, pickle_file_sample_name, model_only_MM, medium_set3, model_only_MM_dict, sample_driven_diet_amounts,\
        sample_sol, post_coop_med, medium_set5, medium_set2, pre_micom_med, pre_micom_med_list, pre_micom_med_df, post_coop_med_df, spent_med_merge, revsigns

########################STEP 3: running human models. 2nd loop.
######SUBTEP 0: file lists and names##########
#have created a new prime metadata file that contains all necessary file names and metadata for this project
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/")
all_metadata= pd.read_csv("Metadata_for_human_and_microbial_communities_Burns_data.csv")
#will need to think about how to split this up into dictionaries for running. BUT, it is most certainly doable. 
#need subdicts for each type of model (tumor/normal)
#addloop#metadata_tumor= all_metadata[all_metadata["Description"]=="tumor"]
#addloop#metadata_healthy= all_metadata[all_metadata["Description"]=="normal"]
#prepare dictionary for looping. NOTE for this example we're running on all six samples, so no need to subset
host_model_dict = dict(zip(sub_metadata.Tissue_Tube_ID, sub_metadata.human_model_file))
#move to correct directory
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_results")

#########create TWO spentmed dicts
healthy_spentmed_file_list = []
healthy_spentmed_file_names = []
for file in os.listdir(): #assuming that listdir contains the correct spentmed files:
    if str("_healthy_microbiome_TD0.2_EUstd_spent_medium.csv") in str(file):
        healthy_spentmed_file = pd.read_csv(file)
        healthy_spentmed_file.drop(columns={"Unnamed: 0"}, inplace=True)
        healthy_spentmed_file_list.append(healthy_spentmed_file)
        healthy_spentmed_file_names.append(str(file))
healthy_spentmed_file_names2 = [x.replace('_healthy_microbiome_TD0.2_EUstd_spent_medium.csv', '') for x in healthy_spentmed_file_names]
healthy_spentmed_dict = dict(zip(healthy_spentmed_file_names2, healthy_spentmed_file_list))
del healthy_spentmed_file_list, healthy_spentmed_file_names, healthy_spentmed_file_names2, healthy_spentmed_file
tumor_spentmed_file_list = []
tumor_spentmed_file_names = []
for file in os.listdir(): #assuming that listdir contains the correct spentmed files:
    if str("_tumor_microbiome_TD0.2_EUstd_spent_medium.csv") in str(file):
        tumor_spentmed_file = pd.read_csv(file)
        tumor_spentmed_file.drop(columns={"Unnamed: 0"}, inplace=True)
        tumor_spentmed_file_list.append(tumor_spentmed_file)
        tumor_spentmed_file_names.append(str(file))
tumor_spentmed_file_names2 = [x.replace('_tumor_microbiome_TD0.2_EUstd_spent_medium.csv', '') for x in tumor_spentmed_file_names]
tumor_spentmed_dict = dict(zip(tumor_spentmed_file_names2, tumor_spentmed_file_list))
del tumor_spentmed_file_list, tumor_spentmed_file_names, tumor_spentmed_file_names2, tumor_spentmed_file
 
#note: we are running each spentmed file TWICE: once on healthy samples, once on tumor samples. 
#split host model dict into two: one for healthy samples, one for tumor samples
healthy_host_model_dict = dict(zip(sub_metadata_normal.Tissue_Tube_ID, sub_metadata_normal.human_model_file))
tumor_host_model_dict = dict(zip(sub_metadata_tumor.Tissue_Tube_ID, sub_metadata_tumor.human_model_file))
#NOW we can match by PBIDs

#create dicts for matching human models to a specific PBID (allows matching to microbiome spent medium files)
host_id_by_PBID_dict_healthy = dict(zip(sub_metadata_normal.Tissue_Tube_ID, sub_metadata_normal.Patient_Blind_ID))
host_id_by_PBID_dict_tumor = dict(zip(sub_metadata_tumor.Tissue_Tube_ID, sub_metadata_tumor.Patient_Blind_ID))
#this is the same as host_to_microbiome_dict?

#create dicts for matching microbiome spent med to a specific PBID; note that we'll be using PBID to find sampleid, so PBIDs should be keys
microbiome_id_by_PBID_dict_healthy = dict(zip(sub_metadata_normal.Patient_Blind_ID, sub_metadata_normal.SampleID))
microbiome_id_by_PBID_dict_tumor = dict(zip(sub_metadata_tumor.Patient_Blind_ID, sub_metadata_tumor.SampleID))

#####FIRST, run healthy spent medium with healthy host models   
for key, value in healthy_host_model_dict.items(): #SUBSET 1: healthy host models
    os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_sample_specific_human_models/")
    human_model = cobra.io.load_matlab_model(value)
    human_model_id = str(key)
    max_growth= human_model.slim_optimize() #run model to get model's MM
    human_model_only_MM =  minimal_medium(human_model, max_growth) #get the model's MM
    human_model_only_MM_dict = human_model_only_MM.to_dict() #creates a dictionary of model's minimal Rxs and their fluxes
    #now, we search for the correct matching spent medium
    healthy_host_PBID = host_id_by_PBID_dict_healthy.get(human_model_id)
    pre_spent_med = healthy_spentmed_dict.get(str("PBID_"+str(healthy_host_PBID))) #SUBSET2: healthy spent medium    
    healthy_spentmed_sampleid = microbiome_id_by_PBID_dict_healthy.get(healthy_host_PBID)#get the sampleid from patient_blind_id
    spent_med = dict(zip(pre_spent_med.premed_name, pre_spent_med.net_postmed_amount)) #make a diet dictionary of metabolites and their amounts
    human_model_exchange_ids = [exchange.id for exchange in human_model.exchanges]  #list comprehension: gets all the export Rxs (exchange ids) from the  model   
    human_to_keep_from_spentmed = {} #initiate a dictionary of metabs/amounts to keep
    for key, value in spent_med.items():
        if key in human_model_exchange_ids:
            human_to_keep_from_spentmed[key]= value
    del key, value
    #get a diet that the sample can grow on
    keys = set(human_to_keep_from_spentmed.keys()).union(human_model_only_MM_dict.keys()) #identify metabolites from diet that can be imported by model
    human_driven_diet_amounts = {k:max(human_to_keep_from_spentmed.get(k,float('-inf')), human_model_only_MM_dict.get(k, float('-inf'))) for k in keys} #merge these dicts
    ###possibly save human_driven_diet_amounts here?
    human_model.medium= human_driven_diet_amounts #set human model medium
    human_sample_sol = human_model.optimize() #run the model
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_results/")
    microbiome_plus_host_name = str("PBID_"+str(healthy_host_PBID)+"_healthy_host_"+str(human_model_id)+"_healthy_microbiome_"\
                                    +str(healthy_spentmed_sampleid)) #this is how we're going to name our files
    fluxes = human_sample_sol.fluxes #pull out fluxes
    #add useful columns: PBID, host_model_id, host_model_description, microbiome-Sampleid, microbiome_description
    fluxes["Patient_blind_ID"] = healthy_host_PBID
    fluxes["host_model_id"] = human_model_id
    fluxes["host_model_description"] = "normal"
    fluxes["microbiome_model_id"] = healthy_spentmed_sampleid
    fluxes["microbiome_model_description"] = "normal"
    fluxes.to_csv("_EUstd_fluxes.csv") #save fluxes... this naming onyl works when we're running with tumor microbiome spent med. 
    human_GR = str(human_sample_sol) #pull out growth rate/biomass
    with open("_EUstd_GRs.csv", "w") as text_file:
        text_file.write(human_GR)
    del text_file
    human_post_growth_med = minimal_medium(human_model, exports=True) #pull out growth medium
    human_post_growth_med["Patient_blind_ID"] = healthy_host_PBID
    human_post_growth_med["host_model_id"] = human_model_id
    human_post_growth_med["host_model_description"] = "normal"
    human_post_growth_med["microbiome_model_id"] = healthy_spentmed_sampleid
    human_post_growth_med["microbiome_model_description"] = "normal"
    human_post_growth_med.to_csv("_uptakes_and_secretions.csv")
    for file in os.listdir():
        src= file 
        if fnmatch.fnmatch(file, "_EUstd_*"): #this should cover all outputs
            dst=str(microbiome_plus_host_name)+file
            os.rename(src,dst)

#####SECOND, run healthy spent medium with tumor host models   
for key, value in tumor_host_model_dict.items(): #SUBSET 1: healthy host models
    os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_sample_specific_human_models/")
    human_model = cobra.io.load_matlab_model(value)
    human_model_id = str(key)
    max_growth= human_model.slim_optimize() #run model to get model's MM
    human_model_only_MM =  minimal_medium(human_model, max_growth) #get the model's MM
    human_model_only_MM_dict = human_model_only_MM.to_dict() #creates a dictionary of model's minimal Rxs and their fluxes
    #now, we search for the correct matching spent medium
    tumor_host_PBID = host_id_by_PBID_dict_tumor.get(human_model_id)
    pre_spent_med = healthy_spentmed_dict.get(str("PBID_"+str(tumor_host_PBID))) #SUBSET2: healthy spent medium    
    healthy_spentmed_sampleid = microbiome_id_by_PBID_dict_healthy.get(tumor_host_PBID)#get the sampleid from patient_blind_id
    spent_med = dict(zip(pre_spent_med.premed_name, pre_spent_med.net_postmed_amount)) #make a diet dictionary of metabolites and their amounts
    human_model_exchange_ids = [exchange.id for exchange in human_model.exchanges]  #list comprehension: gets all the export Rxs (exchange ids) from the  model   
    human_to_keep_from_spentmed = {} #initiate a dictionary of metabs/amounts to keep
    for key, value in spent_med.items():
        if key in human_model_exchange_ids:
            human_to_keep_from_spentmed[key]= value
    del key, value
    #get a diet that the sample can grow on
    keys = set(human_to_keep_from_spentmed.keys()).union(human_model_only_MM_dict.keys()) #identify metabolites from diet that can be imported by model
    human_driven_diet_amounts = {k:max(human_to_keep_from_spentmed.get(k,float('-inf')), human_model_only_MM_dict.get(k, float('-inf'))) for k in keys} #merge these dicts
    ###possibly save human_driven_diet_amounts here?
    human_model.medium= human_driven_diet_amounts #set human model medium
    human_sample_sol = human_model.optimize() #run the model
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_results/")
    microbiome_plus_host_name = str("PBID_"+str(tumor_host_PBID)+"_tumor_host_"+str(human_model_id)+"_healthy_microbiome_"\
                                    +str(healthy_spentmed_sampleid)) #this is how we're going to name our files
    fluxes = human_sample_sol.fluxes #pull out fluxes
    #add useful columns: PBID, host_model_id, host_model_description, microbiome-Sampleid, microbiome_description
    fluxes.to_csv("_EUstd_fluxes.csv") #save fluxes... this naming onyl works when we're running with tumor microbiome spent med. 
    human_GR = str(human_sample_sol) #pull out growth rate/biomass
    with open("_EUstd_GRs.csv", "w") as text_file:
        text_file.write(human_GR)
    del text_file
    human_post_growth_med = minimal_medium(human_model, exports=True) #pull out growth medium
    human_post_growth_med.to_csv("_uptakes_and_secretions.csv")
    for file in os.listdir():
        src= file 
        if fnmatch.fnmatch(file, "_EUstd_*"): #this should cover all outputs
            dst=str(microbiome_plus_host_name)+file
            os.rename(src,dst)

#####THIRD, run tumor spent medium with healthy host models   
for key, value in healthy_host_model_dict.items(): #SUBSET 1: healthy host models
    os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_sample_specific_human_models/")
    human_model = cobra.io.load_matlab_model(value)
    human_model_id = str(key)
    max_growth= human_model.slim_optimize() #run model to get model's MM
    human_model_only_MM =  minimal_medium(human_model, max_growth) #get the model's MM
    human_model_only_MM_dict = human_model_only_MM.to_dict() #creates a dictionary of model's minimal Rxs and their fluxes
    #now, we search for the correct matching spent medium
    healthy_host_PBID = host_id_by_PBID_dict_healthy.get(human_model_id)
    pre_spent_med = tumor_spentmed_dict.get(str("PBID_"+str(healthy_host_PBID))) #SUBSET2: healthy spent medium    
    tumor_spentmed_sampleid = microbiome_id_by_PBID_dict_tumor.get(healthy_host_PBID)#get the sampleid from patient_blind_id
    spent_med = dict(zip(pre_spent_med.premed_name, pre_spent_med.net_postmed_amount)) #make a diet dictionary of metabolites and their amounts
    human_model_exchange_ids = [exchange.id for exchange in human_model.exchanges]  #list comprehension: gets all the export Rxs (exchange ids) from the  model   
    human_to_keep_from_spentmed = {} #initiate a dictionary of metabs/amounts to keep
    for key, value in spent_med.items():
        if key in human_model_exchange_ids:
            human_to_keep_from_spentmed[key]= value
    del key, value
    #get a diet that the sample can grow on
    keys = set(human_to_keep_from_spentmed.keys()).union(human_model_only_MM_dict.keys()) #identify metabolites from diet that can be imported by model
    human_driven_diet_amounts = {k:max(human_to_keep_from_spentmed.get(k,float('-inf')), human_model_only_MM_dict.get(k, float('-inf'))) for k in keys} #merge these dicts
    ###possibly save human_driven_diet_amounts here?
    human_model.medium= human_driven_diet_amounts #set human model medium
    human_sample_sol = human_model.optimize() #run the model
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_results/")
    microbiome_plus_host_name = str("PBID_"+str(healthy_host_PBID)+"_healthy_host_"+str(human_model_id)+"_tumor_microbiome_"\
                                    +str(tumor_spentmed_sampleid)) #this is how we're going to name our files
    fluxes = human_sample_sol.fluxes #pull out fluxes
    #add useful columns: PBID, host_model_id, host_model_description, microbiome-Sampleid, microbiome_description
    fluxes.to_csv("_EUstd_fluxes.csv") #save fluxes... this naming onyl works when we're running with tumor microbiome spent med. 
    human_GR = str(human_sample_sol) #pull out growth rate/biomass
    with open("_EUstd_GRs.csv", "w") as text_file:
        text_file.write(human_GR)
    del text_file
    human_post_growth_med = minimal_medium(human_model, exports=True) #pull out growth medium
    human_post_growth_med.to_csv("_uptakes_and_secretions.csv")
    for file in os.listdir():
        src= file 
        if fnmatch.fnmatch(file, "_EUstd_*"): #this should cover all outputs
            dst=str(microbiome_plus_host_name)+file
            os.rename(src,dst)

#####FOURTH, run tumor spent medium with tumor host models   
for key, value in tumor_host_model_dict.items(): #SUBSET 1: healthy host models
    os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_sample_specific_human_models/")
    human_model = cobra.io.load_matlab_model(value)
    human_model_id = str(key)
    max_growth= human_model.slim_optimize() #run model to get model's MM
    human_model_only_MM =  minimal_medium(human_model, max_growth) #get the model's MM
    human_model_only_MM_dict = human_model_only_MM.to_dict() #creates a dictionary of model's minimal Rxs and their fluxes
    #now, we search for the correct matching spent medium
    tumor_host_PBID = host_id_by_PBID_dict_tumor.get(human_model_id)
    pre_spent_med = tumor_spentmed_dict.get(str("PBID_"+str(tumor_host_PBID))) #SUBSET2: healthy spent medium    
    tumor_spentmed_sampleid = microbiome_id_by_PBID_dict_tumor.get(tumor_host_PBID)#get the sampleid from patient_blind_id
    spent_med = dict(zip(pre_spent_med.premed_name, pre_spent_med.net_postmed_amount)) #make a diet dictionary of metabolites and their amounts
    human_model_exchange_ids = [exchange.id for exchange in human_model.exchanges]  #list comprehension: gets all the export Rxs (exchange ids) from the  model   
    human_to_keep_from_spentmed = {} #initiate a dictionary of metabs/amounts to keep
    for key, value in spent_med.items():
        if key in human_model_exchange_ids:
            human_to_keep_from_spentmed[key]= value
    del key, value
    #get a diet that the sample can grow on
    keys = set(human_to_keep_from_spentmed.keys()).union(human_model_only_MM_dict.keys()) #identify metabolites from diet that can be imported by model
    human_driven_diet_amounts = {k:max(human_to_keep_from_spentmed.get(k,float('-inf')), human_model_only_MM_dict.get(k, float('-inf'))) for k in keys} #merge these dicts
    ###possibly save human_driven_diet_amounts here?
    human_model.medium= human_driven_diet_amounts #set human model medium
    human_sample_sol = human_model.optimize() #run the model
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_diet_results/")
    microbiome_plus_host_name = str("PBID_"+str(tumor_host_PBID)+"_tumor_host_"+str(human_model_id)+"_tumor_microbiome_"\
                                    +str(tumor_spentmed_sampleid)) #this is how we're going to name our files
    fluxes = human_sample_sol.fluxes #pull out fluxes
    #add useful columns: PBID, host_model_id, host_model_description, microbiome-Sampleid, microbiome_description
    fluxes.to_csv("_EUstd_fluxes.csv") #save fluxes... this naming onyl works when we're running with tumor microbiome spent med. 
    human_GR = str(human_sample_sol) #pull out growth rate/biomass
    with open("_EUstd_GRs.csv", "w") as text_file:
        text_file.write(human_GR)
    del text_file
    human_post_growth_med = minimal_medium(human_model, exports=True) #pull out growth medium
    human_post_growth_med.to_csv("_uptakes_and_secretions.csv")
    for file in os.listdir():
        src= file 
        if fnmatch.fnmatch(file, "_EUstd_*"): #this should cover all outputs
            dst=str(microbiome_plus_host_name)+file
            os.rename(src,dst)




























