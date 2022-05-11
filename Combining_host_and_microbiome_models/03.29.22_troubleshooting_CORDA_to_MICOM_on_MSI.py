#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 17:14:32 2022

@author: adamo010
"""

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
from itertools import chain, groupby
import cobra
import micom
from micom import load_pickle
from micom.media import minimal_medium

###############STEP 0: file lists and names################ #we are ONLY running normal microbiome and tumor human samples here today.
def main():
    #script = sys.argv[0]
    microbiome_model = "Sample_05_community.pickle" #these files are ONLY tumor human files. see Burns2015_tumor_pickle_file_list.txt
    microbiome_model_sample_prename = str("Sample_05_community.pickle")
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/")
    all_metadata= pd.read_csv("Metadata_for_human_and_microbial_communities_Burns_data_V2.csv")
    metadata_subset= all_metadata[all_metadata["microbiome_model_file"]==microbiome_model_sample_prename] #CHANGE WITH EACH script
    #host_model_dict = dict(zip(metadata_subset.host_sampleid, metadata_subset.human_model_file))
    #host_model_names_dict = dict(zip(metadata_subset.host_sampleid, metadata_subset.host_description))
    relevant_PBID = str(metadata_subset.iloc[1,1])   
    #microbiome_model_dict = dict(zip(metadata_subset.microbiome_sampleid, metadata_subset.microbiome_model_file))
    #host_type_dict = dict(zip(metadata_subset.microbiome_sampelid, metadata_subset. ))
    #os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Merging_Burns_MICOM_and_CORDA/CORDA_spent_med_output_round2/")
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Combining_host_and_microbiome_models/EUstd_host_spent_med_localgen/")
    metadata_subset["pbid_host_desc"] = metadata_subset["Patient_Blind_ID"].astype(str)+ "_"+ metadata_subset["host_description"].astype(str)
    host_sampleinfo = dict(zip(metadata_subset.pbid_host_desc, metadata_subset.host_sampleid))
    spentmed_file_list = []
    spentmed_file_names = [] #CHANGE _normal/tumor_hostmodel_EU.. .csv below when using tumor/normal host
    for file in os.listdir():
    	if str(relevant_PBID) in str(file):
    		spentmed_file = pd.read_csv(file)
    		spentmed_file.drop(columns={"Unnamed: 0"}, inplace=True)
            #spentmed_file_sampleid = host_sampleinfo.get(relevant_PBID( ) 
    		spentmed_file_list.append(spentmed_file)
    		spentmed_file_names.append(str(file))
    spentmed_file_names = [x.replace('_hostmodel_EUstd_spent_medium.csv', '') for x in spentmed_file_names] #CHANGE ME when using tumor/normal host
    spentmed_file_names = [x.replace('PBID_', '') for x in spentmed_file_names]
    spentmed_dict= dict(zip(spentmed_file_names, spentmed_file_list))
    #microbiome_id_by_PBID_dict = dict(zip(metadata_subset.microbiome_sampleid, metadata_subset.Patient_Blind_ID))
    #host_id_by_PBID_dict = dict(zip(metadata_subset.host_sampleid, metadata_subset.Patient_Blind_ID)) #edited
    #os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Patient_specific_microbiomes/")
    microbiome_model = load_pickle(microbiome_model) #import our pickle file
    microbiome_model_id = microbiome_model_sample_prename.replace("_community.pickle","") #get model from file name
    microbiome_model_only_MM =  minimal_medium(microbiome_model, 1.0) #get model_only_MM
    microbiome_model_only_MM_dict = microbiome_model_only_MM.to_dict() #creates a dictionary of model's minimal Rxs and their fluxes  
    #NEXT STEP: figure out how to extract host model/spent med file 
    for key, value in spentmed_dict.items():
        host_spent_med = value
        host_proxy_sampleid = str(key)
        host_real_sampleid  = host_sampleinfo.get(host_proxy_sampleid) #returns the actual sampleid, vs PBID_desc
        host_type = host_proxy_sampleid.split("_")[1] #HAHA! assign host type to variable. Less changing by hand!
        #microbiome_PBID = microbiome_id_by_PBID_dict.get(microbiome_model_id) #now, we search for the correct matching spent medium
        #pre_spent_med = spentmed_dict.get(str("PBID_"+str(relevant_PBID)))   
        #spentmed_sampleid = host_id_by_PBID_dict.get(microbiome_PBID) #get the host sampleid from patient_Blind_id
        host_spent_med_dict = dict(zip(host_spent_med.premed_name, host_spent_med.net_postmed_amount)) #make a diet dictionary of metabolites and their amounts
        host_spent_med_dict = {k.replace("[e]", "_m"): v for k, v in host_spent_med.items()} #need to edit human model metab names to gel with microbiomes       
        microbiome_model_exchange_ids = [exchange.id for exchange in microbiome_model.exchanges]  #list comprehension: gets all the export Rxs (exchange ids) from the  model   
        microbiome_to_keep_from_spentmed = {} #initiate a dictionary of metabs/amounts to keep
        for key, value in host_spent_med_dict.items():
            if key in microbiome_model_exchange_ids:
                microbiome_to_keep_from_spentmed[key]= value
        #del key, value
        keys = set(microbiome_to_keep_from_spentmed.keys()).union(microbiome_model_only_MM_dict.keys()) #identify metabolites from diet that can be imported by model
        microbiome_driven_diet_amounts = {k:max(microbiome_to_keep_from_spentmed.get(k,float('-inf')), microbiome_model_only_MM_dict.get(k, float('-inf'))) for k in keys} #merge these dicts
        microbiome_model.medium= microbiome_driven_diet_amounts #set microbiome model medium
        microbiome_sample_sol = microbiome_model.cooperative_tradeoff(fraction=0.2, fluxes=True, pfba=False) #run FBA with TD=0.2 FOR THIS FILE ONLY turn fluxes off
        os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Merging_Burns_MICOM_and_CORDA/MICOM_spent_med_output_round2/")
        host_plus_microbiome_name = "PBID_"+str(relevant_PBID)+"_"+str(microbiome_model_id)+"_tumor_microbiome_"+"to_"+str(host_real_sampleid)+\
                                    "_"+str(host_type)+"_host" #this is how we're going to name our files; #CHANGE WITH EACH script
        fluxes = microbiome_sample_sol.fluxes #pull out fluxes
        fluxes.rename(columns={"fluxes": "flux_amount"}, inplace=True)
        fluxes.to_csv(host_plus_microbiome_name+"_EUstd_fluxes.csv") #save fluxes.
        rates = microbiome_sample_sol.members.growth_rate.drop("medium")
        rates= rates.to_frame()
        rates["Patient_Blind_ID"] = relevant_PBID
        rates["microbiome_model_id"] = microbiome_model_id
        rates["microbiome_model_description"] = "tumor" #CHANGE WITH EACH script
        rates["host_model_id"] = host_real_sampleid
        rates["host_model_description"] = host_type
        rates.to_csv(host_plus_microbiome_name+"_EUstd_indiv_spp_GRs.csv") #save fluxes.
        microbiome_GR = str(microbiome_sample_sol) #pull out growth rate/biomass
        with open(host_plus_microbiome_name+"_EUstd_GRs.csv", "w") as text_file:
            text_file.write(microbiome_GR)
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Merging_Burns_MICOM_and_CORDA/")        	        	
    return

main()
