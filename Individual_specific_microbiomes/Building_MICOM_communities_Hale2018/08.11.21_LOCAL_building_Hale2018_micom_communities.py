#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 13:56:36 2021

@author: adamo010
"""
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

#I know I'm going to regret this, but here goes- building Hale pickle models locally. 

os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Building_MICOM_communities_Hale2018/")

def micom_for_sample(sample_id_file):
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Building_MICOM_communities_Hale2018/Hale_OTU_specific_abundance_tables/")
    abundance_dataframe = pd.read_csv(str(sample_id_file))
    abundance_file_sample_name = str(sample_id_file).rstrip("_abundances_and_models.csv")
    abundance_dataframe.rename(columns = {'OTU_ID': 'id', abundance_file_sample_name: 'abundance', 'otu_file_name': 'file'}, inplace=True)
    abundance_dataframe['id']= abundance_dataframe['id'].astype(str)
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/Hale_2018_model_building_and_combining/Hale_2018_otu_specific_models/")
    sample_community = Community(abundance_dataframe, solver='gurobi')
    os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Building_MICOM_communities_Hale2018/testo/")
    sample_community.to_pickle("_community.pickle")
    for file in os.listdir():
    	src= file
    	if fnmatch.fnmatch(file, "_community.pickle"):
    		dst = str(abundance_file_sample_name)+file
    		os.rename(src, dst)
    del file, src, dst, sample_community	     
    return 

micom_for_sample("s007F13xB1A01_abundances_and_models.csv")     


os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Building_MICOM_communities_Hale2018/Hale_OTU_specific_abundance_tables/")
abundance_dataframe = pd.read_csv("s007F13xB1A01_abundances_and_models.csv")








       
