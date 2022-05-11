#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 15:35:31 2020

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

#pick up where we left off from MICOM_healthy_vs_CRC_trial.py: saved the solutions as pickle files

#another thing I would like to do is see if I can compare the growth media of different
#communities. Do inputs/outputs differ? Goal: make a graph of the differences 

from micom import load_pickle
healthycomm = load_pickle("healthy_5spp_community.pickle")
CRCcomm = load_pickle("CRC_5spp_community.pickle")
#these are the built communities. Use this to run FBA, cooperative FBA, etc. 

healthyFBA = healthycomm.optimize()
CRCFBA = CRCcomm.optimize()

healthyFBA.members
CRCFBA.members
    
from micom.media import minimal_medium

healthycomm_MM = minimal_medium(healthycomm, 0.8)    #note that minimizing components is unfeasibly slow (i.e. won't run)
CRCcomm_MM = minimal_medium(CRCcomm, 0.8)

#ok, making minimal media is fun. Is it useful?
healthycomm_MM2 = minimal_medium(healthycomm, 0.8, exports=True)
CRCcomm_MM2 = minimal_medium(CRCcomm, 0.8, exports=True)
#now it is: exports are negative, imports are positive. This is whole-community import/export at 80% max community growth rate
#to do tom: try different max community growth rates, and compare healthy and diseased.

possible_GRs = np.arange(0.1,1.01,0.1)
possible_GRs = np.around(possible_GRs, decimals=1)      #this takes care of the weird multiple decimal values that np gives for 0.3 and 0.7

##########HEALTHY#################
#create the initial data frame to be added to later
healthycomm_MM_GR = minimal_medium(healthycomm, 1.0, exports=True)
healthycomm_MM_GR= healthycomm_MM_GR.reset_index()
healthycomm_MM_GR.rename(columns={"index": "Exchange_Rx", 0: ("Flux"+"_"+"1.0")}, inplace=True)    #rename columns appropriately

for elem in possible_GRs:
    healthycomm_MM_temp = minimal_medium(healthycomm, elem, exports=True)   #run the model with the appropriate GR value and save it as a temp.pd
    healthycomm_MM_temp= healthycomm_MM_temp.reset_index()  #reset model.pd column index
    healthycomm_MM_temp.rename(columns={"index": "Exchange_Rx", 0: (str(elem)+"_"+"Flux")}, inplace=True)    #rename columns appropriately
    healthycomm_MM_GR = pd.merge(left = healthycomm_MM_temp, right= healthycomm_MM_GR, how= "left", left_on="Exchange_Rx", right_on="Exchange_Rx")
    #merge model.pd with a running pandas series containing exchange Rxs 
healthycomm_MM_GR= healthycomm_MM_GR.drop("Flux_1.0", axis=1)       #clean up the column crom the initial data frame
#very cool!

#now, we must edit the exchange_Rx names
#######from Fuso_diets_making_graphs
metab_key = pd.read_csv("Metabolites_key_full.tsv", sep='\t')    #read in the metabolite names key
for col in metab_key.columns:
    print(col)
#there are a lot of extra columns here: ony really want the first 2
metab_key_slimmed = metab_key.loc[:, 'abbreviation':'fullName']    #use : to select full column
metab_key_slimmed.head()    #cool
metab_key_slimmed = metab_key_slimmed.rename(columns={"abbreviation": "Abbreviation",
                                                     "fullName": "FullName"})
#code below takes all the metabolite long names and puts them in a dictionary with short names    
metab_key_slimmed2 = metab_key_slimmed.copy()
metab_key_slimmed2.head()
type(metab_key_slimmed2)
metabolite_list_dict = metab_key_slimmed2.set_index("Abbreviation").to_dict()["FullName"]
print(metabolite_list_dict)
#double check
for key in metabolite_list_dict.keys():
    print(key)
    
#from healthycomm_MM_GR, need to remove EX_ and _m from front and back of items in Exchange_Rx column, respectively 
healthycomm_MM_GR["Exchange_Rx"] = healthycomm_MM_GR["Exchange_Rx"].str.replace("EX_", "")
healthycomm_MM_GR
healthycomm_MM_GR["Exchange_Rx"] = healthycomm_MM_GR["Exchange_Rx"].str.replace("_m", "")
healthycomm_MM_GR
#excellent! Now, to replace short forms in Exchange_Rx with long forms

healthycomm_MM_GR_fullnames = copy.deepcopy(healthycomm_MM_GR)
#want to do a vlookup kind of thing: replace short form names in Exchange_Rx column with long form names in metabolite_list_dict
healthycomm_MM_GR_fullnames['Exchange_Rx']= healthycomm_MM_GR_fullnames['Exchange_Rx'].map(metabolite_list_dict)
#well, that's a lot easier than changing column names

print(healthycomm_MM_GR_fullnames.dtypes)   #we are good on the variable type front.

#graphing:
SMALL_SIZE = 15
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

healthycomm_MM_GR_fullnames.plot.bar(figsize=(10,6), rot=0) 
plt.xlabel("Metabolites")
plt.ylabel("Flux value, mmol/gDCW/h")
#ok so this graph is terrible and also doesn't include metabolite list names
#make two graphs: one for import and one for export.
#basically, scan through each row/metabolite name and check for negative values
#negative values get added to a new dataframe called 'healthy_efflux'
#positive values get added to a new dataframe called 'healthy_uptake'

#first, check that there are no discordant rows (i.e. where import switches to export or vice versa across different community GRs)
#you know what? I'm not in the mood for this shit. Just look at the dataframe. It'll change colour if the numbers are discordant. Here, they're fine.

uptake_export = []
for row in healthycomm_MM_GR_fullnames["1.0_Flux"]:
    if row < 0:
        uptake_export.append("export")
    else:
        uptake_export.append("uptake")

healthycomm_MM_GR_fullnames['Uptake/Export']= uptake_export        
#create a list, append export or uptake to that list based on whether each value in column 1.0_flux
#in the pandas series healthycomm_MM_GR_fullnames, then add that column. -ve values in 1.0_flux 
#should match 'export' entries and +ve values should match 'import'

healthy_export = healthycomm_MM_GR_fullnames[healthycomm_MM_GR_fullnames["Uptake/Export"]=="export"]
print(healthy_export)
healthy_uptake = healthycomm_MM_GR_fullnames[healthycomm_MM_GR_fullnames["Uptake/Export"]=="uptake"]
print(healthy_uptake)
#create new pandas series healthy_export and healthy_uptake based on values from Uptake/Export column

healthy_export.set_index("Exchange_Rx", inplace=True)
healthy_uptake.set_index("Exchange_Rx", inplace=True)
#have to reset index in each or stupid matplotlib will just number the x-axis instead of providing the x-axis labels

healthy_export.plot.bar(figsize=(10,6), rot=0) 
plt.xlabel("Metabolite", labelpad = 15)
plt.xticks(rotation=45)
plt.ylabel("Flux value, mmol/gDCW/h", labelpad = 15)
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.title("MICOM- predicted healthy 5-species community exports", size=BIGGER_SIZE)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5)) #this moves legend to outside of the graph

healthy_uptake.plot.bar(figsize=(10,6), rot=0) 
plt.xlabel("Metabolite", labelpad = 15)
plt.xticks(rotation=90)
plt.ylabel("Flux value, mmol/gDCW/h", labelpad = 15)
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.title("MICOM- predicted healthy 5-species community uptakes", size=BIGGER_SIZE)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5)) #this moves legend to outside of the graph

#well, these aren't ideal, but they're fine.

##########CRC#################
#create the initial data frame to be added to later
CRCcomm_MM_GR = minimal_medium(CRCcomm, 1.0, exports=True)
CRCcomm_MM_GR= CRCcomm_MM_GR.reset_index()
CRCcomm_MM_GR.rename(columns={"index": "Exchange_Rx", 0: ("Flux"+"_"+"1.0")}, inplace=True)    #rename columns appropriately

for elem in possible_GRs:
    CRCcomm_MM_temp = minimal_medium(CRCcomm, elem, exports=True)   #run the model with the appropriate GR value and save it as a temp.pd
    CRCcomm_MM_temp= CRCcomm_MM_temp.reset_index()  #reset model.pd column index
    CRCcomm_MM_temp.rename(columns={"index": "Exchange_Rx", 0: (str(elem)+"_"+"Flux")}, inplace=True)    #rename columns appropriately
    CRCcomm_MM_GR = pd.merge(left = CRCcomm_MM_temp, right= CRCcomm_MM_GR, how= "left", left_on="Exchange_Rx", right_on="Exchange_Rx")
    #merge model.pd with a running pandas series containing exchange Rxs 
CRCcomm_MM_GR= CRCcomm_MM_GR.drop("Flux_1.0", axis=1)       #clean up the column crom the initial data frame
#very cool!

#now, we must edit the exchange_Rx names
#######from Fuso_diets_making_graphs
metab_key = pd.read_csv("Metabolites_key_full.tsv", sep='\t')    #read in the metabolite names key
for col in metab_key.columns:
    print(col)
#there are a lot of extra columns here: ony really want the first 2
metab_key_slimmed = metab_key.loc[:, 'abbreviation':'fullName']    #use : to select full column
metab_key_slimmed.head()    #cool
metab_key_slimmed = metab_key_slimmed.rename(columns={"abbreviation": "Abbreviation",
                                                     "fullName": "FullName"})
#code below takes all the metabolite long names and puts them in a dictionary with short names    
metab_key_slimmed2 = metab_key_slimmed.copy()
metab_key_slimmed2.head()
type(metab_key_slimmed2)
metabolite_list_dict = metab_key_slimmed2.set_index("Abbreviation").to_dict()["FullName"]
print(metabolite_list_dict)
#double check
for key in metabolite_list_dict.keys():
    print(key)
    
#from CRCcomm_MM_GR, need to remove EX_ and _m from front and back of items in Exchange_Rx column, respectively 
CRCcomm_MM_GR["Exchange_Rx"] = CRCcomm_MM_GR["Exchange_Rx"].str.replace("EX_", "")
CRCcomm_MM_GR
CRCcomm_MM_GR["Exchange_Rx"] = CRCcomm_MM_GR["Exchange_Rx"].str.replace("_m", "")
CRCcomm_MM_GR
#excellent! Now, to replace short forms in Exchange_Rx with long forms

CRCcomm_MM_GR_fullnames = copy.deepcopy(CRCcomm_MM_GR)
#want to do a vlookup kind of thing: replace short form names in Exchange_Rx column with long form names in metabolite_list_dict
CRCcomm_MM_GR_fullnames['Exchange_Rx']= CRCcomm_MM_GR_fullnames['Exchange_Rx'].map(metabolite_list_dict)
#well, that's a lot easier than changing column names

print(CRCcomm_MM_GR_fullnames.dtypes)   #we are good on the variable type front.

#graphing:
SMALL_SIZE = 15
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

CRCcomm_MM_GR_fullnames.plot.bar(figsize=(10,6), rot=0) 
plt.xlabel("Metabolites")
plt.ylabel("Flux value, mmol/gDCW/h")
#ok so this graph is terrible and also doesn't include metabolite list names
#make two graphs: one for import and one for export.
#basically, scan through each row/metabolite name and check for negative values
#negative values get added to a new dataframe called 'CRC_efflux'
#positive values get added to a new dataframe called 'CRC_uptake'

#first, check that there are no discordant rows (i.e. where import switches to export or vice versa across different community GRs)
#you know what? I'm not in the mood for this shit. Just look at the dataframe. It'll change colour if the numbers are discordant. Here, they're fine.

uptake_export = []
for row in CRCcomm_MM_GR_fullnames["1.0_Flux"]:
    if row < 0:
        uptake_export.append("export")
    else:
        uptake_export.append("uptake")

CRCcomm_MM_GR_fullnames['Uptake/Export']= uptake_export        
#create a list, append export or uptake to that list based on whether each value in column 1.0_flux
#in the pandas series CRCcomm_MM_GR_fullnames, then add that column. -ve values in 1.0_flux 
#should match 'export' entries and +ve values should match 'import'

CRC_export = CRCcomm_MM_GR_fullnames[CRCcomm_MM_GR_fullnames["Uptake/Export"]=="export"]
print(CRC_export)
CRC_export=CRC_export.replace("Protein-Linked Serine Or Threonine Residue (O-Glycosylation Site)", "Serine/Threonine protein")
#added
for row in CRC_export["Exchange_Rx"]:
    print (row)
    
CRC_uptake = CRCcomm_MM_GR_fullnames[CRCcomm_MM_GR_fullnames["Uptake/Export"]=="uptake"]
print(CRC_uptake)
#create new pandas series CRC_export and CRC_uptake based on values from Uptake/Export column

CRC_export.set_index("Exchange_Rx", inplace=True)
CRC_uptake.set_index("Exchange_Rx", inplace=True)
#have to reset index in each or stupid matplotlib will just number the x-axis instead of providing the x-axis labels

CRC_export.plot.bar(figsize=(10,6), rot=0) 
plt.xlabel("Metabolite", labelpad = 15)
plt.xticks(rotation=45)
plt.ylabel("Flux value, mmol/gDCW/h", labelpad = 15)
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.title("MICOM- predicted CRC 5-species community exports", size=BIGGER_SIZE)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5)) #this moves legend to outside of the graph

CRC_uptake.plot.bar(figsize=(10,6), rot=0) 
plt.xlabel("Metabolite", labelpad = 15)
plt.xticks(rotation=90)
plt.ylabel("Flux value, mmol/gDCW/h", labelpad = 15)
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.title("MICOM- predicted CRC 5-species community uptakes", size=BIGGER_SIZE)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5)) #this moves legend to outside of the graph

#well, these aren't ideal, but they're fine.
    






















