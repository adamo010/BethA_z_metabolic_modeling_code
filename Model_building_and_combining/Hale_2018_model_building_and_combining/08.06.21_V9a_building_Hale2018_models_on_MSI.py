import sys
import micom
from micom.util import (
    load_model,
    join_models,
    add_var_from_expression,
    adjust_solver_config,
    clean_ids,
)
import gurobipy
import pickle
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
import ast
from ast import literal_eval

#the goal of this code is to combine all the models together where multiple models have been assigned to a single OTU

#we will operate out of the directory Combining_Hale_models

def main():
    script = sys.argv[0]
    agora_file = sys.argv[1]
    all_agora_merged_condensed2= pd.read_pickle(agora_file)
    all_agora_merged_condensed2 = all_agora_merged_condensed2.drop(columns=["OTU_ID"])
    agora_model_file_dict = dict(zip(all_agora_merged_condensed2.OTU_proxy_ID,all_agora_merged_condensed2.all_agora_model_file_names.apply(lambda x: literal_eval(str(x)))))    
    agora_dict_A = dict(list(agora_model_file_dict.items())[len(agora_model_file_dict)//2:])
    agora_dict_B = dict(list(agora_model_file_dict.items())[:len(agora_model_file_dict)//2])
    agora_dict_1 = dict(list(agora_dict_A.items())[len(agora_dict_A)//2:])
    agora_dict_2 = dict(list(agora_dict_A.items())[:len(agora_dict_A)//2])
    agora_dict_3 = dict(list(agora_dict_B.items())[len(agora_dict_B)//2:])
    agora_dict_4 = dict(list(agora_dict_B.items())[:len(agora_dict_B)//2])
    del agora_dict_A, agora_dict_B  
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/AGORA_models/")
    for key, value in agora_dict_1.items():
        filename = str(key)+"_model.xml"
        otu_model = micom.util.join_models(value, str(key))
        os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/OTU_models")
        cobra.io.write_sbml_model(otu_model, filename)
        del otu_model
        os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/AGORA_models/")
    del key, value, filename
    for key, value in agora_dict_2.items():
        filename = str(key)+"_model.xml"
        otu_model = micom.util.join_models(value, str(key))
        os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/OTU_models")
        cobra.io.write_sbml_model(otu_model, filename)
        del otu_model
        os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/AGORA_models/")
    del key, value, filename
    for key, value in agora_dict_3.items():
        filename = str(key)+"_model.xml"
        otu_model = micom.util.join_models(value, str(key))
        os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/OTU_models")
        cobra.io.write_sbml_model(otu_model, filename)
        del otu_model
        os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/AGORA_models/")
    del key, value, filename
    for key, value in agora_dict_4.items():
        filename = str(key)+"_model.xml"
        otu_model = micom.util.join_models(value, str(key))
        os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/OTU_models")
        cobra.io.write_sbml_model(otu_model, filename)
        del otu_model
        os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/AGORA_models/")
    del key, value, filename
    return

main()    




