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

###############STEP 0: file lists and names################
def main():
    script = sys.argv[0]
    pickle_file = sys.argv[1] #these files are ONLY tumor microbiome files. 
    pickle_file_sample_prename = str(sys.argv[1])
    all_metadata= pd.read_csv("Metadata_for_human_and_microbial_communities_Burns_data.csv")
    sub_metadata_tumor = all_metadata[all_metadata["Description"]=="tumor"]
    #sub_metadata_normal = all_metadata[all_metadata["Description"]=="normal"]
    pickle_file_sample_name = pickle_file_sample_prename.replace("_community.pickle", "")
    EUstandard = {} #editing with every iteration, so need to reset.
    with open("EUstandard_fluxes.tsv", "r") as file:
    	next(file)
    	for line in file:                                   
        	linestuff = line.strip().split()                    
        	EUstandard[linestuff[0]] = float(linestuff[1])  #this creates a dictionary called 'EUstandard' containing the metabolite (key) and its amount (value) from the supplied diet file
    EUstandard_edited = {k.replace("[e]", "_m"): v for k, v in EUstandard.items()}
    EUstandard_edited2 = {k.replace("(e)", "_m"): v for k, v in EUstandard_edited.items()}
    del file, line, linestuff, EUstandard, EUstandard_edited
    EU_diet = EUstandard_edited2.items()
    EU_diet_list = list(EU_diet)
    medium_set1 = pd.DataFrame(EU_diet_list)
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Patient_specific_microbiomes")
    #normal_microbiome_model_dict = dict(zip(sub_metadata_normal.SampleID, sub_metadata_normal.microbiome_model_file))
    tumor_microbiome_model_dict = dict(zip(sub_metadata_tumor.SampleID, sub_metadata_tumor.microbiome_model_file))
    sample_to_PBID_dict = dict(zip(all_metadata.SampleID, all_metadata.Patient_Blind_ID)) #key is sampleid (unique), value is PBID (not unique)
    imported_pickle_model = load_pickle(pickle_file) #import our pickle file
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
    del useable_EU_diet, useable_EU_diet_listc    keys = set(to_keep_from_EU_diet.keys()).union(model_only_MM_dict.keys()) #create a union of model MM and model-useable EUstd diet keys
    sample_driven_diet_amounts = {k:max(to_keep_from_EU_diet.get(k,float('-inf')), model_only_MM_dict.get(k, float('-inf'))) for k in keys} #merge dicts
    microbiome_premed = sample_driven_diet_amounts.items() #save medium information for later
    microbiome_premed_list = list(microbiome_premed)
    medium_set4 = pd.DataFrame(microbiome_premed_list)
    del microbiome_premed, microbiome_premed_list, keys
    ########Save any important medium output files
    imported_pickle_model.medium= sample_driven_diet_amounts #set sample_driven_diet_amounts to model's growth medium
    sample_sol = imported_pickle_model.cooperative_tradeoff(fraction=0.3, fluxes=False, pfba=False) #run FBA with TD=0.2 FOR THIS FILE ONLY turn fluxes off
    post_coop_med = minimal_medium(imported_pickle_model, sample_sol.growth_rate, exports=True)#now, extract spent medium from this.
    #medium_set5 = copy.deepcopy(post_coop_med).to_frame() #save medium information for later
    #medium_set5.reset_index(inplace=True) #clean up the dataframe
    pre_micom_med = sample_driven_diet_amounts.items() #convert sample_driven_diet_amounts to a dataframe.
    pre_micom_med_list = list(pre_micom_med)
    pre_micom_med_df = pd.DataFrame(pre_micom_med_list)
    pre_micom_med_df.rename(columns={0: "premed_name", 1: "premed_amount"}, inplace=True) #then, rename columns etc.
    post_coop_med_df = post_coop_med.to_frame() #convert post_coop_med (which is a series) to dataframe
    post_coop_med_df.reset_index(inplace=True)
    post_coop_med_df.rename(columns={"index": "postmed_name", 0: "postmed_amount"}, inplace=True)
    spent_med_merge= pd.merge(left=pre_micom_med_df, right=post_coop_med_df, how="left", left_on="premed_name", right_on="postmed_name") 
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
    medium_set7a = pd.merge(left = medium_set3, right= medium_set1, how='left', left_on="index", right_on=0, indicator=True)
    medium_set7 = medium_set7a[medium_set7a['_merge'] == "left_only"] 
    medium_set7.drop(columns={"0_y", 1, "_merge"}, inplace=True)
    del medium_set7a
    medium_set8a = pd.merge(left = medium_set1, right= medium_set2, how='left', left_on=0, right_on=0, indicator=True)
    medium_set8 = medium_set8a[medium_set8a['_merge'] == "left_only"] 
    medium_set8.drop(columns={"1_y", "_merge"}, inplace=True)
    del medium_set8a
    #now, we want to get spent medium, which is medium set 6-7+8
    medium_set9a = pd.merge(left=medium_set6, right= medium_set7, how="left", left_on="premed_name", right_on="index", indicator=True)
    medium_set9b = medium_set9a[medium_set9a["_merge"] != "both"]
    medium_set9b.drop(columns={"premed_amount", "index", "0_x", "_merge"}, inplace=True)
    #great. Next step- add metabolites from set 8 (part of EUstd but not useable by the model. )
    medium_set8.rename(columns={0: "premed_name", "1_x": "net_postmed_amount"}, inplace=True) #stack EUstd diets to merge vertically
    medium_set9= pd.concat([medium_set9b, medium_set8], axis=0) #now, merge for the final one!
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Merging_Burns_MICOM_and_CORDA/MICOM_spent_med_output")
    medium_set9.to_csv("_TD0.2_EUstd_spent_medium.csv")  #save this medium
    pbid = sample_to_PBID_dict.get(pickle_file_sample_name) #ectract PBID (value) from sample_to_PBID_dict
    unique_file_name=  "PBID_"+ str(pbid)+"_tumor_microbiome"
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_TD0.2_EUstd_spent_medium.csv"):
            dst=str(unique_file_name)+file
            os.rename(src,dst)
    return

main()    
