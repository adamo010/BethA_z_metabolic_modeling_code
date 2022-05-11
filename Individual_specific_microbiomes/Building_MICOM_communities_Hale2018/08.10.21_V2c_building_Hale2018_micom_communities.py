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
    abundance_file = sys.argv[1]
    abundance_file_sample_prename = str(sys.argv[1])
    abundance_file_sample_name = abundance_file_sample_prename.replace("_abundances_and_models.csv", "")
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/OTU_abundances_and_models/")	
    abundance_dataframe = pd.read_csv(abundance_file)
    abundance_dataframe.rename(columns = {'OTU_ID': 'id', abundance_file_sample_name: 'abundance', 'otu_file_name': 'file'}, inplace=True)
    abundance_dataframe['id']= abundance_dataframe['id'].astype(str)
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/OTU_models") #directory where otu models live
    sample_community = Community(abundance_dataframe, solver='gurobi')
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/Hale2018_pickle_files/")
    sample_community.to_pickle("_community.pickle")
    for file in os.listdir():
    	src= file
    	if fnmatch.fnmatch(file, "_community.pickle"):
    		dst = str(abundance_file_sample_name)+file
    		os.rename(src, dst)
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/OTU_models/")		
    del key, value, file, src, dst, sample_community	     
    return

main()    


