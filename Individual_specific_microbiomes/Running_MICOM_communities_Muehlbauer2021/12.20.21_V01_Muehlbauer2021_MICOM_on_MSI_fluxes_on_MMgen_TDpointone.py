import os
import copy
import pandas as pd
import corda
from corda import reaction_confidence
from cobra.io import load_matlab_model
from cobra.io import read_sbml_model
from corda import CORDA    
import gurobipy
import scipy
import numpy as np
import cobra
import csv
import subprocess
import re
from itertools import chain, groupby
import fnmatch
import sys #for specifying input files externally

#we start in adamo010/CORDA_Burns_data/

def main():
    script = sys.argv[0]
    GPR_file = sys.argv[1]
    GPR_file_sample_prename = str(sys.argv[1])
    GPR_file_sample_name = GPR_file_sample_prename.replace("_GPRs_binned_V3.csv", "")
    recon3 = load_matlab_model("/panfs/roc/groups/7/blekhman/adamo010/CORDA_Burns_data/Recon3DModel_301.mat")
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/CORDA_Burns_data/CORDA_GPR_files/")
    #first, import and clean up GPR files
    imported_GPR_file = pd.read_csv(GPR_file)
    imported_GPR_file.drop(columns={"Unnamed: 0", "edited_Rx_rule"}, inplace=True)
    abb_Rx_name = []
    for row in imported_GPR_file["Rx_name"]:
    	splitlist = row.split(r": ") #remove everything after :
    	newrow = splitlist[0]
    	abb_Rx_name.append(newrow)
    del row, splitlist, newrow
    imported_GPR_file["abb_Rx_name"] = abb_Rx_name
    imported_GPR_file.drop(columns=['Rx_name'], inplace=True)
    gene_conf = dict(zip(imported_GPR_file.abb_Rx_name, imported_GPR_file.solved_Rx_rule))
	#second, CORDA-ize the model
    opt = CORDA(recon3, gene_conf)
    opt.build() #RATE LIMITING STEP
    opt.cobra_model() #convert to cobra model
    opt.model.objective="biomass_maintenance"
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/CORDA_Burns_data/Completed_CORDA_models/")
    cobra.io.save_matlab_model(opt.model, "_V3_model.mat")
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_V3_model.mat"):
            dst=str(GPR_file_sample_name)+file
            os.rename(src,dst)
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/CORDA_Burns_data/")
    return

main()