#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 14:56:37 2020

@author: adamo010
"""
#my previous approach missed a tonne of taxa, for some reason. Need to figure out why and fix.
#so, at least with Akkermansia, the problem appears to be based on matching whole phylogenies vs species names.
#eg. in AGORA, Akkermansia muciniphilia is in the Akkermansiaceae family; in Burns data, it's in Verrucomicrobiaceae
#how do I fix this without making a mess? Can't just match on genus names. 
#would matching based on, e.g., genus, vs genus_level_sort, be more helpful?


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

#I don't want to re-invent the wheel, so I'm going to use the output from Building_microbiome_communities_03.25.20.py

#re-do agora model list
agora = pd.read_csv("AGORA_model_phylogeny_table.csv")
adapted_agora = copy.deepcopy(agora)  #copy to keep everything clean

#I could create the taxonomy_level_sort stuff again, but I don't think that's helpful
#instead, I'm only going to create two new columns: genus/species, and genus2/species
#then, I will create the same two new columns in the otu table and match on those.

adapted_agora["genus_species"] = "g__" + adapted_agora['genus'] + ";  s__" + adapted_agora['species']
adapted_agora["genus_species2"] = "g__" + adapted_agora['genus2'] + ";  s__" + adapted_agora['species']

adapted_agora['genus_species'] = adapted_agora['genus_species'].str.strip()
adapted_agora['genus_species2'] = adapted_agora['genus_species2'].str.strip()

#now, I will create g-s models dictionaries where all models with the same genus_species id are compiled in a list
unique_agora_genus_spp = adapted_agora.genus_species.unique()
unique_agora_genus_spp = pd.DataFrame(unique_agora_genus_spp)
unique_agora_genus_spp['model_list']=np.empty((len(unique_agora_genus_spp), 0)).tolist()
unique_agora_genus_spp.reset_index()
unique_agora_genus_spp = unique_agora_genus_spp.rename(columns={0: "genus_species"})

unique_agora_genus_spp2 = adapted_agora.genus_species2.unique()
unique_agora_genus_spp2 = pd.DataFrame(unique_agora_genus_spp2)
unique_agora_genus_spp2['model_list']=np.empty((len(unique_agora_genus_spp2), 0)).tolist()
unique_agora_genus_spp2.reset_index()
unique_agora_genus_spp2 = unique_agora_genus_spp2.rename(columns={0: "genus_species2"})

#now, have two dataframes ready to input lists of matching models into. 

genus_spp_models_dict = {}
for row in unique_agora_genus_spp['genus_species']:
    gs_name_key = str(row)
    temp_df =adapted_agora[adapted_agora["genus_species"].str.contains(str(row))]
    templist = []
    for row in temp_df["model_ID"]:
        templist.append(str(row))
    del temp_df
    genus_spp_models_dict[gs_name_key] = templist
    del templist
    del gs_name_key

genus_spp2_models_dict = {}
for row in unique_agora_genus_spp2['genus_species2']:
    gs2_name_key = str(row)
    temp_df =adapted_agora[adapted_agora["genus_species2"].str.contains(str(row))]
    templist = []
    for row in temp_df["model_ID"]:
        templist.append(str(row))
    del temp_df
    genus_spp2_models_dict[gs2_name_key] = templist
    del templist
    del gs2_name_key

#cool. Let's import the model-free otu tables
modelless_otus_V1 = pd.read_csv("03.26.20_OTUs_without_models_V1.csv")
modelless_otus_V2 = pd.read_csv("03.26.20_OTUs_without_models_V2.csv")

#add a new column to match on the newly created dictionaries
modelless_otus_V1["genus_species"] = modelless_otus_V1['genus'] + "; " + modelless_otus_V1['species']
modelless_otus_V2["genus_species2"] = modelless_otus_V2['genus'] + "; " + modelless_otus_V2['species']

modelless_otus_V1["genus_species"] = modelless_otus_V1["genus_species"].str.strip()
modelless_otus_V2["genus_species2"] = modelless_otus_V2["genus_species2"].str.strip()

modelless_otus_V1['genus_species_model_list'] = modelless_otus_V1['genus_species'].map(genus_spp_models_dict)        
modelless_otus_V2['genus_species2_model_list'] = modelless_otus_V2['genus_species2'].map(genus_spp2_models_dict)        

#couple things: 
#1) in this case, at least, cutting the white space mattered
#2) this added like three model-assigned OTUs. Wahoo

#pull out the assigned otus
otus_with_gs_models_V1 = modelless_otus_V1[modelless_otus_V1['genus_species_model_list'].notnull()]   
otus_with_gs_models_V2 = modelless_otus_V2[modelless_otus_V2['genus_species2_model_list'].notnull()]   
    
#import previously assigned model list
prev_assigned_V1 = pd.read_csv("03.26.20_OTUs_with_models_V1.csv")
prev_assigned_V2 = pd.read_csv("03.26.20_OTUs_with_models_V2.csv")

#concatenate (2step process)
V1 = [prev_assigned_V1, otus_with_gs_models_V1]
all_V1 = pd.concat(V1)
V2 = [prev_assigned_V2, otus_with_gs_models_V2]
all_V2 = pd.concat(V2)

all_V1.to_csv("03.27.20_OTUs_with_models_V1.csv")
all_V2.to_csv("03.27.20_OTUs_with_models_V2.csv")

#pull out the aotus that still don't have assigned models
otus_wo_gs_models_V1 = modelless_otus_V1[modelless_otus_V1['genus_species_model_list'].isnull()]   
otus_wo_gs_models_V2 = modelless_otus_V2[modelless_otus_V2['genus_species2_model_list'].isnull()]   
    
otus_wo_gs_models_V1.to_csv("03.27.20_OTUs_without_models_V1.csv")
otus_wo_gs_models_V2.to_csv("03.27.20_OTUs_without_models_V2.csv")

#great. So now we're back to finding new models.


















