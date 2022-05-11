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

#we will operate out of the directory Combining_Muehlbauer_models

def main():
    script = sys.argv[0]
    carveme_file = sys.argv[1]
    all_carveme_merged_condensed2= pd.read_pickle(carveme_file)
    #all_agora_merged_condensed2 = all_agora_merged_condensed2.drop(columns=["OTU_ID"])
    carveme_model_file_dict = dict(zip(all_carveme_merged_condensed2.OTU_ID,all_carveme_merged_condensed2.all_carveme_matched_models.apply(lambda x: literal_eval(str(x)))))    
    os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Muehlbauer_models/CARVEME_models/")
    for key, value in carveme_model_file_dict.items():
        filename = str(key)+"_model.xml"
        otu_model = micom.util.join_models(value, str(key))
        os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Muehlbauer_models/OTU_models")
        cobra.io.write_sbml_model(otu_model, filename)
        del otu_model
        os.chdir("/panfs/roc/groups/7/blekhman/adamo010/Combining_Muehlbauer_models/CARVEME_models/")
    return

main()    





