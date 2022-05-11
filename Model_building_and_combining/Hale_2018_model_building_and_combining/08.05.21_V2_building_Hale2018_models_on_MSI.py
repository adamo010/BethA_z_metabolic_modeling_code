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
    all_agora_merged_condensed2 = pd.read_pickle('AGORA_models_that_need_OTU_builds.pickle')
    all_agora_merged_condensed2 = all_agora_merged_condensed2.drop(columns=["OTU_ID"])
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/AGORA_models/")
    agora_model_file_dict = dict(zip(all_agora_merged_condensed2.OTU_proxy_ID, \
                                 all_agora_merged_condensed2.all_agora_model_file_names.apply(lambda x: literal_eval(str(x)))))    
    for key, value in agora_model_file_dict.items():
        filename = str(key)+"_model.xml"
        otu_model = micom.util.join_models(value, str(key))
        os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/OTU_models")
        cobra.io.write_sbml_model(otu_model, filename)
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/")
    all_mambo_merged_condensed2 = pd.read_pickle('MAMBO_models_that_need_OTU_builds.pickle')
    all_mambo_merged_condensed2= all_mambo_merged_condensed2.drop(columns=["OTU_ID"])
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/MAMBO_models/")
    mambo_model_file_dict = dict(zip(all_mambo_merged_condensed2.OTU_proxy_ID, \
                                 all_mambo_merged_condensed2.all_mambo_model_file_names.apply(lambda x: literal_eval(str(x)))))    
    for key, value in mambo_model_file_dict.items():
        filename = str(key)+"_model.xml"
        otu_model = micom.util.join_models(value, str(key))
        os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/OTU_models")
        cobra.io.write_sbml_model(otu_model, filename)
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/")
    carveme2 = pd.read_pickle('CarveMe_models_that_need_OTU_builds.pickle')
    carveme2= carveme2.drop(columns=["OTU_ID"])
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/CarveMe_models/")
    carveme_model_file_dict = dict(zip(carveme2.OTU_proxy_ID, \
                                 carveme2.all_mambo_model_file_names.apply(lambda x: literal_eval(str(x)))))    
    for key, value in carveme_model_file_dict.items():
        filename = str(key)+"_model.xml"
        otu_model = micom.util.join_models(value, str(key))
        os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/OTU_models")
        cobra.io.write_sbml_model(otu_model, filename)
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Hale_models/")



