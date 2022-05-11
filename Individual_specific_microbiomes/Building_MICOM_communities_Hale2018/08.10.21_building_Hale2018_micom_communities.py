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
from micom import Community    

def main():
    script = sys.argv[0]
    abundance_files = sys.argv[1]
    my_file = open(abundance_files, "r")
    abundance_files = my_file.readlines()
    #step 1: create a list of all the files with each sample's OTU table
    sample_OTU_table_files = []
    for file in glob.glob("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/OTU_abundances_and_models/*abundances_and_models.csv"):
    	if file in abundance_files:
    		sample_OTU_abundance_files.append(file)
    del file
    #step 2: create a list of sample names
    sample_OTU_table_samplenames = []
    for elem in sample_OTU_table_files:
    	elem = elem.replace("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/OTU_abundances_and_models/","")
    	elem = elem.replace("_abundances_and_models.csv","")
    	sample_OTU_table_samplenames.append(elem)
    del elem
    #step 3: merge file names and sample names into a dictionary.
    sample_dict = {sample_OTU_table_samplenames[i]: sample_OTU_table_files[i] for i in range(len(sample_OTU_table_samplenames))}
    #step 4: import and edit all dataframes to match requirements for Community function in micom
    edited_sample_dict = {}
    for key, value in sample_dict.items():
    	dataframe = pd.read_csv(value)
    	dataframe.rename(columns = {'OTU_ID': 'id', key: 'abundance', 'otu_file_name': 'file'}, inplace=True)
    	dataframe['id']= dataframe['id'].astype(str)
    	edited_sample_dict.update( {key : dataframe})
    del key, value
    #step 5: run micom and save pickle files
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/OTU_models") #directory where otu models live
    finished_communities = {}
    for key, value in edited_sample_dict.items():
    	sample_community = Community(value, solver='gurobi')
    	os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/Hale2018_pickle_files/")
    	sample_community.to_pickle("_community.pickle")
    	for file in os.listdir():
    		src= file
    		if fnmatch.fnmatch(file, "_community.pickle"):
    			dst = str(key)+file
    			os.rename(src, dst)
    	os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/OTU_models/")		
    del key, value, file, src, dst	     
    return

main()    


