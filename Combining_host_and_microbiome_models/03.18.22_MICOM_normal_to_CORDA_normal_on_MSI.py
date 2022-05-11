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
from cobra.io import load_matlab_model
from cobra.medium import minimal_medium
import corda
from corda import CORDA

###############STEP 0: file lists and names################ 
#we are ONLY running normal microbiome and normal human samples here today.
def main():
    script = sys.argv[0]
    human_model = sys.argv[1] #these files are ONLY normal human files. 
    human_model_sample_prename = str(sys.argv[1])
    all_metadata= pd.read_csv("Metadata_for_human_and_microbial_communities_Burns_data_V2.csv")
    metadata_subset= all_metadata[all_metadata["Pair_type"]=="normal_normal"] #CHANGE WITH EACH script
    host_model_dict = dict(zip(metadata_subset.host_sampleid, metadata_subset.human_model_file))
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Merging_Burns_MICOM_and_CORDA/MICOM_spent_med_output/")
    spentmed_file_list = []
    spentmed_file_names = []
    for file in os.listdir():
    	if str("_normal_microbiome_TD0.2_EUstd_spent_medium.csv") in str(file): #CHANGE ME when using tumor/normal microbiome
        	spentmed_file = pd.read_csv(file)
        	spentmed_file.drop(columns={"Unnamed: 0"}, inplace=True)
        	spentmed_file_list.append(spentmed_file)
        	spentmed_file_names.append(str(file))
    spentmed_file_names = [x.replace('_normal_microbiome_TD0.2_EUstd_spent_medium.csv', '') for x in spentmed_file_names] #CHANGE ME when using tumor/normal microbiome
    spentmed_dict = dict(zip(spentmed_file_names, spentmed_file_list))
    host_id_by_PBID_dict = dict(zip(metadata_subset.host_sampleid, metadata_subset.Patient_Blind_ID))
    microbiome_id_by_PBID_dict = dict(zip(metadata_subset.Patient_Blind_ID, metadata_subset.microbiome_sampleid)) #edited
    for key, value in host_model_dict.items():
        os.chdir("/panfs/roc/groups/7/blekhman/adamo010/CORDA_Burns_data/Completed_CORDA_models/")
        human_model = cobra.io.load_matlab_model(value)
        human_model_id = str(key)
        max_growth= human_model.slim_optimize() #run model to get model's MM
        human_model_only_MM =  minimal_medium(human_model, max_growth) #get the model's MM
        human_model_only_MM_dict = human_model_only_MM.to_dict() #creates a dictionary of model's minimal Rxs and their fluxes
        host_PBID = host_id_by_PBID_dict.get(human_model_id) #now, we search for the correct matching spent medium
        pre_spent_med = spentmed_dict.get(str("PBID_"+str(host_PBID)))   
        spentmed_sampleid = microbiome_id_by_PBID_dict.get(host_PBID)#get the sampleid from patient_blind_id
        spent_med = dict(zip(pre_spent_med.premed_name, pre_spent_med.net_postmed_amount)) #make a diet dictionary of metabolites and their amounts
        human_model_exchange_ids = [exchange.id for exchange in human_model.exchanges]  #list comprehension: gets all the export Rxs (exchange ids) from the  model   
        human_to_keep_from_spentmed = {} #initiate a dictionary of metabs/amounts to keep
        for key, value in spent_med.items():
            if key in human_model_exchange_ids:
                human_to_keep_from_spentmed[key]= value
        #del key, value
        keys = set(human_to_keep_from_spentmed.keys()).union(human_model_only_MM_dict.keys()) #identify metabolites from diet that can be imported by model
        human_driven_diet_amounts = {k:max(human_to_keep_from_spentmed.get(k,float('-inf')), human_model_only_MM_dict.get(k, float('-inf'))) for k in keys} #merge these dicts
        human_model.medium= human_driven_diet_amounts #set human model medium
        human_sample_sol = human_model.optimize() #run the model
        os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Merging_Burns_MICOM_and_CORDA/CORDA_spent_med_output")
        microbiome_plus_host_name = str("PBID_"+str(host_PBID)+"_normal_host_"+str(human_model_id)+"_normal_microbiome_"\
                                    +str(spentmed_sampleid)) #this is how we're going to name our files; #CHANGE WITH EACH script
        fluxes = human_sample_sol.fluxes #pull out fluxes
        fluxes = fluxes.to_frame()
        fluxes.rename(columns={"fluxes": "flux_amount"}, inplace=True)
        fluxes["Patient_blind_ID"] = host_PBID
        fluxes["host_model_id"] = human_model_id
        fluxes["host_model_description"] = "normal" #CHANGE WITH EACH script
        fluxes["microbiome_model_id"] = spentmed_sampleid
        fluxes["microbiome_model_description"] = "normal" #CHANGE WITH EACH script
        fluxes.to_csv(microbiome_plus_host_name+"_EUstd_fluxes.csv") #save fluxes.
        human_GR = str(human_sample_sol) #pull out growth rate/biomass
        with open(microbiome_plus_host_name+"_EUstd_GRs.csv", "w") as text_file:
            text_file.write(human_GR)
        #del text_file
        human_post_growth_med = minimal_medium(human_model, exports=True) #pull out growth medium
        human_post_growth_med = human_post_growth_med.to_frame() #convert series to dataframe
        human_post_growth_med.rename(columns={0: "metab_amount"}, inplace=True)
        human_post_growth_med["Patient_blind_ID"] = host_PBID
        human_post_growth_med["host_model_id"] = human_model_id
        human_post_growth_med["host_model_description"] = "normal" #CHANGE WITH EACH script
        human_post_growth_med["microbiome_model_id"] = spentmed_sampleid
        human_post_growth_med["microbiome_model_description"] = "normal" #CHANGE WITH EACH script
        human_post_growth_med.to_csv(microbiome_plus_host_name+"_EUstd_uptakes_and_secretions.csv")
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Merging_Burns_MICOM_and_CORDA/")        	        	
    return

main()


        
    


    


