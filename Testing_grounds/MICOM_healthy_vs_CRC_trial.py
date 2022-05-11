#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 16:16:14 2020

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

#OK, let's try MICOM! First thing we need to do is create a pandas series with 
#in /Users/adamo010/Documents/microbiome_FBA

###step 1: create the community
#at least two columns: ID (which specifies the taxon ID) and file (which specifies where the model is)

healthyspp_files = {"Fprau": "/Users/adamo010/Documents/microbiome_FBA/Fprau.json", 
              "Bbif": "/Users/adamo010/Documents/microbiome_FBA/Bbif.json",
              "Ccat": "/Users/adamo010/Documents/microbiome_FBA/Ccat.json",
              "Lacid": "/Users/adamo010/Documents/microbiome_FBA/Lacid.json",
              "Rinul": "/Users/adamo010/Documents/microbiome_FBA/Rinul.json"}

#side note: use cobrapy to get information on the number of reactions, metabolites, and genes

Fprau_model = cobra.io.load_json_model("Fprau.json")
Bbif_model = cobra.io.load_json_model("Bbif.json")
Ccat_model = cobra.io.load_json_model("Ccat.json")
Lacid_model = cobra.io.load_json_model("Lacid.json")
Rinul_model = cobra.io.load_json_model("Rinul.json")

healthy_model_list = [Fprau_model, Bbif_model, Ccat_model, Lacid_model, Rinul_model]
healthy_model_dict = {"Fprau": str(Fprau_model), "Bbif": str(Bbif_model),
                      "Ccat": str(Ccat_model), "Lacid": str(Lacid_model), "Rinul": str(Rinul_model)}

healthy_num_reactions = {}
healthy_num_metabolites = {}
healthy_num_genes = {}

for model in healthy_model_list:
    #print(str(model))
    #print(len(model.reactions))
    healthy_num_reactions.update({str(model): len(model.reactions)})
    #print(len(model.metabolites))
    healthy_num_metabolites.update({str(model): len(model.metabolites)})
    #print(len(model.genes)) 
    healthy_num_genes.update({str(model): len(model.genes)})

####at this point, need to merge several dictionaries:
    #healthyspp2 (contains id and file)
    #healthy_model_dict (contains id and full model name)
    #healthy_num_reactions2, healthy_num_metabolites2, and healthy_num_genes2: contain #s and full model names
# so I think step 1 is to merge healthy_model_dict and healthy_num_dicts so #s are associated with id rather than full model name
#do this in a pandas series

healthyspp = pd.Series(healthy_model_dict)
healthyspp = healthyspp.reset_index()
healthyspp.columns = ["id", "model_name"]    
healthyspp['reactions']= healthyspp['model_name'].map(healthy_num_reactions)
healthyspp['metabolites']= healthyspp['model_name'].map(healthy_num_metabolites)
healthyspp['genes']= healthyspp['model_name'].map(healthy_num_genes)
healthyspp['file'] = healthyspp['id'].map(healthyspp_files)

###step 2: community growth rates and fluxes
from micom import Community #this function creates a community from a pandas series of models
healthycomm_FBA = Community(healthyspp, solver='gurobi')

from micom import load_pickle
healthycomm_FBA.to_pickle("healthy_5spp_community.pickle")

#now, we have run MICOM FBA for the models in healthycomm and stored the results in healthycomm_FBA

print(healthycomm_FBA.objective.expression)
#this prints the objective expression for MICOM; by default it is maximising community growth rate

%time healthysol1 = healthycomm_FBA.optimize() #prints out the community model solution along with how long it took to arrive there
healthysol1.members  #prints out the community members' abundance, growth rates, Rxs, and metabolites

%time healthysol2 = healthycomm_FBA.optimize(fluxes=True) #prints out the community model solution along with how long it took to arrive there, plus fluxes
healthysol2.fluxes

#####this is cool, I can work with this. 

###step 3: COOPERATIVE TRADEOFF!!!

%time healthyspp_CT = healthycomm_FBA.cooperative_tradeoff(fraction=1.0)
#got a stupid "GLPK only supports linear objectives" error. Why isn't MICOM using Gurobi???
#naturally, now that I restarted the kernel, it magically works again.
healthyspp_CT

#set some growth rate limits
healthyspp_CT_min1 = healthycomm_FBA.cooperative_tradeoff(min_growth=0.1)  # single value
healthyspp_CT_min1
healthyspp_CT_min2 = healthycomm_FBA.cooperative_tradeoff(min_growth=[0.1, 0.2, 0.3, 0.4, 0.5])
healthyspp_CT_min2

#test the impact of the tradeoff parameter on your solution
#make sure numpy is engaged
%time healthyspp_CT_tradeoff = healthycomm_FBA.cooperative_tradeoff(fraction=np.arange(0.1, 1.01, 0.1))

#tried to save the model as a pickle file and it didn't work. probably b/c i'm saving a variety of solutions, not just one
healthyspp_CT_tradeoff.to_pickle("healthy_coop_tradeoff.pickle")
from micom import load_pickle
healthyspp_CT_tradeoff2 = load_pickle("healthy_coop_tradeoff.pickle")
#THIS IS THE COOLEST FEATURE but PSA takes a long time to run

#at its core, waht is healthyspp_CT_tradeoff? A pandas series containing 2 columns. One has the tradeoff values, the other has the solutions
#Presumably, 'the usual pandas means' referred to in the tutorial means looking at each solution



#######CLEANUP TIME!!!############
try1 = copy.deepcopy(healthyspp_CT_tradeoff)
try1["solution"] = try1["solution"].astype(str)
try1["solution"] = try1["solution"].str.replace("<CommunitySolution ", "")  #remove first part of string. 
try1["solution"] = try1["solution"].str.slice(stop=-16) #remove the last 16 characters from this string
try1["solution"] = try1["solution"].str.strip() #clear white space from both
#try1["tradeoff"] = try1["tradeoff"].str.strip() #clear white space from both
try1["solution"] = try1["solution"].astype(float)   #convert values to float for graphing purposes
#try1.round({"tradeoff": 1}) #this will hopefully eliminate the need for the next line.
try1.tradeoff = [1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]   #Not all the tradeoff values are rounded- this didn't work
try1.rename(columns={"solution": "Community", 'tradeoff': "Tradeoff"}, inplace=True)    #rename columns appropriately

healthy_tradeoff_cleaned = copy.deepcopy(try1)

#ok, now I have a fun solution. What the hell was I planning on doing with it?

#"the solution can then be inspected by the usual pandas means" that everyone knows duh
#look at individual spp growth rates

healthy_tradeoff_cleaned

healthyspp_CT_tradeoff

healthyspp_rates =healthyspp_CT_tradeoff.solution.apply(lambda x: x.members.growth_rate)
#what this does is goes through the solution column in healthyspp_CT_tradeoff and iterates through each column.
#for each value of solution, apply the lambda function (x represents the solution): get the members.growth rate
#another way to do this:
for row in healthyspp_CT_tradeoff['solution']:
    print(row)
    #print(row.members)      #prints each community members' abundance, growth rate, Rxs, and metabolites
    print(row.members.growth_rate)  #prints each community members' growth rate only 
    #print(row.members.abundance)
    #print(row.members.reactions)
    #print(row.members.metabolites)


healthyspp_rates
del healthyspp_rates['medium']   #remove medium column b/c it's all 'na'

#I think this would be the most fun thing to graph
#also, in a lucky turn of events, 'rates' is already a pandas dataframe- no more dumbass printing to dataframe for me!
#but, I do need to add a column to reflect community growth rate
************************************
healthyspp_rates['Tradeoff'] = [1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1] #I'd really like some way to associate compartment with the tradeoff value

###now, merge rates on shared 'tradeoff' column
healthyspp_tradeoff_all = pd.merge(left= healthyspp_rates, right= healthy_tradeoff_cleaned, how="left",left_on="Tradeoff", right_on="Tradeoff")
healthyspp_tradeoff_all

#now, reorder the columns
healthyspp_tradeoff_all = healthyspp_tradeoff_all[["Tradeoff", "Bbif", "Ccat", "Fprau", "Lacid", "Rinul", "Community"]]

#NEAT

#now, I would like to graph these. X-axis is tradeoff, y-axis is GR
#note: everything needs to be converted to floats for this to actually work
print(healthyspp_tradeoff_all.dtypes)       #aha, just Community is an object(?)

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

#this works- may need to tweak sizes etc.
plt.plot('Tradeoff', 'Community', data=healthyspp_tradeoff_all, marker='o', color='black', markersize = 6, linewidth = 3)
plt.plot('Tradeoff', 'Fprau', data=healthyspp_tradeoff_all, marker='o', color='orange', markersize = 6, linewidth = 3)
plt.plot('Tradeoff', 'Bbif', data=healthyspp_tradeoff_all, marker='o', color='limegreen', markersize = 6, linewidth = 3)
plt.plot('Tradeoff', 'Ccat', data=healthyspp_tradeoff_all, marker='o', color='blue', markersize = 6, linewidth = 3)
plt.plot('Tradeoff', 'Lacid', data=healthyspp_tradeoff_all, marker='o', color='darkviolet', markersize = 6, linewidth = 3)
plt.plot('Tradeoff', 'Rinul', data=healthyspp_tradeoff_all, marker='o', color='crimson', markersize = 6, linewidth = 3)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5)) #this moves legend to outside of the graph
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.xlabel("Tradeoff value", labelpad = 15)
plt.ylabel("Growth rate", labelpad = 15)
plt.title("MICOM- predicted individual and community growth \n rates for five species of healthy-associated microbes", size=BIGGER_SIZE)
#once things are ready, do this:
plt.savefig("/Users/adamo010/Documents/microbiome_FBA/healthy_indiv_plus_comm_GRs_01.28.20.png", dpi=300, bbox_inches="tight", pad_inches=0.1)


###############CRC models#############
###step 1: create the community
#at least two columns: ID (which specifies the taxon ID) and file (which specifies where the model is)

CRCspp_files = {"Amuc": "/Users/adamo010/Documents/microbiome_FBA/Amuc.json", 
              "Bfrag": "/Users/adamo010/Documents/microbiome_FBA/Bfrag.json",
              "Ecoli": "/Users/adamo010/Documents/microbiome_FBA/Ecoli.json",
              "Fnuc": "/Users/adamo010/Documents/microbiome_FBA/Fnuc.json",
              "Pana": "/Users/adamo010/Documents/microbiome_FBA/Pana.json"}

#side note: use cobrapy to get information on the number of reactions, metabolites, and genes

Amuc_model = cobra.io.load_json_model("Amuc.json")
Bfrag_model = cobra.io.load_json_model("Bfrag.json")
Ecoli_model = cobra.io.load_json_model("Ecoli.json")
Fnuc_model = cobra.io.load_json_model("Fnuc.json")
Pana_model = cobra.io.load_json_model("Pana.json")

CRC_model_list = [Amuc_model, Bfrag_model, Ecoli_model, Fnuc_model, Pana_model]
CRC_model_dict = {"Amuc": str(Amuc_model), "Bfrag": str(Bfrag_model),
                      "Ecoli": str(Ecoli_model), "Fnuc": str(Fnuc_model), "Pana": str(Pana_model)}

CRC_num_reactions = {}
CRC_num_metabolites = {}
CRC_num_genes = {}

for model in CRC_model_list:
    #print(str(model))
    #print(len(model.reactions))
    CRC_num_reactions.update({str(model): len(model.reactions)})
    #print(len(model.metabolites))
    CRC_num_metabolites.update({str(model): len(model.metabolites)})
    #print(len(model.genes)) 
    CRC_num_genes.update({str(model): len(model.genes)})

CRCspp = pd.Series(CRC_model_dict)
CRCspp = CRCspp.reset_index()
CRCspp.columns = ["id", "model_name"]    
CRCspp['reactions']= CRCspp['model_name'].map(CRC_num_reactions)
CRCspp['metabolites']= CRCspp['model_name'].map(CRC_num_metabolites)
CRCspp['genes']= CRCspp['model_name'].map(CRC_num_genes)
CRCspp['file'] = CRCspp['id'].map(CRCspp_files)

###step 2: community growth rates and fluxes
from micom import Community #this function creates a community from a pandas series of models
CRCcomm_FBA = Community(CRCspp, solver='gurobi')
#now, we have run MICOM FBA for the models in CRCcomm and stored the results in CRCcomm_FBA

from micom import load_pickle
CRCcomm_FBA.to_pickle("CRC_5spp_community.pickle")


print(CRCcomm_FBA.objective.expression)
#this prints the objective expression for MICOM; by default it is maximising community growth rate

%time CRCsol1 = CRCcomm_FBA.optimize() #prints out the community model solution along with how long it took to arrive there
CRCsol1.members  #prints out the community members' abundance, growth rates, Rxs, and metabolites

%time CRCsol2 = CRCcomm_FBA.optimize(fluxes=True) #prints out the community model solution along with how long it took to arrive there, plus fluxes
CRCsol2.fluxes

###step 3: COOPERATIVE TRADEOFF!!!
%time CRCspp_CT = CRCcomm_FBA.cooperative_tradeoff(fraction=1.0)
#got a stupid "GLPK only supports linear objectives" error. Why isn't MICOM using Gurobi???
#naturally, now that I restarted the kernel, it magically works again.
CRCspp_CT

#set some growth rate limits
CRCspp_CT_min1 = CRCcomm_FBA.cooperative_tradeoff(min_growth=0.1)  # single value
CRCspp_CT_min1
CRCspp_CT_min2 = CRCcomm_FBA.cooperative_tradeoff(min_growth=[0.1, 0.2, 0.3, 0.4, 0.5])
CRCspp_CT_min2

#test the impact of the tradeoff parameter on your solution
#make sure numpy is engaged
%time CRCspp_CT_tradeoff = CRCcomm_FBA.cooperative_tradeoff(fraction=np.arange(0.1, 1.01, 0.1))
CRCspp_CT_tradeoff
#THIS IS THE COOLEST FEATURE but PSA takes a long itme to run

#######CLEANUP TIME!!!############
try1 = copy.deepcopy(CRCspp_CT_tradeoff)
try1["solution"] = try1["solution"].astype(str)
try1["solution"] = try1["solution"].str.replace("<CommunitySolution ", "")  #remove first part of string. 
try1["solution"] = try1["solution"].str.slice(stop=-16) #remove the last 16 characters from this string
try1["solution"] = try1["solution"].str.strip() #clear white space from both
#try1["tradeoff"] = try1["tradeoff"].str.strip() #clear white space from both
try1["solution"] = try1["solution"].astype(float)   #convert values to float for graphing purposes
#try1.round({"tradeoff": 1}) #this will hopefully eliminate the need for the next line.
try1.tradeoff = [1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]   #Not all the tradeoff values are rounded
try1.rename(columns={"solution": "Community", 'tradeoff': "Tradeoff"}, inplace=True)    #rename columns appropriately

CRC_tradeoff_cleaned = copy.deepcopy(try1)

#ok, now I have a fun solution. What the hell was I planning on doing with it?

#"the solution can then be inspected by the usual pandas means" that everyone knows duh
#look at individual spp growth rates
CRCspp_rates =CRCspp_CT_tradeoff.solution.apply(lambda x: x.members.growth_rate)
CRCspp_rates
del CRCspp_rates['medium']   #remove medium column b/c it's all 'na'

#I think this would be the most fun thing to graph
#also, in a lucky turn of events, 'rates' is already a pandas dataframe- no more dumbass printing to dataframe for me!
#but, I do need to add a column to reflect community growth rate
************************************
CRCspp_rates['Tradeoff'] = [1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1] #I'd really like some way to associate compartment with the tradeoff value

###now, merge rates on shared 'tradeoff' column
CRCspp_tradeoff_all = pd.merge(left= CRCspp_rates, right= CRC_tradeoff_cleaned, how="left",left_on="Tradeoff", right_on="Tradeoff")
CRCspp_tradeoff_all

#now, reorder the columns
CRCspp_tradeoff_all = CRCspp_tradeoff_all[["Tradeoff", "Amuc", "Bfrag", "Ecoli", "Fnuc", "Pana", "Community"]]

#NEAT

#now, I would like to graph these. X-axis is tradeoff, y-axis is GR
#note: everything needs to be converted to floats for this to actually work
print(CRCspp_tradeoff_all.dtypes)       #aha, just Community is an object(?)


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

#this works- may need to tweak sizes etc.
plt.plot('Tradeoff', 'Community', data=CRCspp_tradeoff_all, marker='o', color='black', markersize = 6, linewidth = 3)
plt.plot('Tradeoff', 'Amuc', data=CRCspp_tradeoff_all, marker='o', color='orange', markersize = 6, linewidth = 3)
plt.plot('Tradeoff', 'Bfrag', data=CRCspp_tradeoff_all, marker='o', color='limegreen', markersize = 6, linewidth = 3)
plt.plot('Tradeoff', 'Ecoli', data=CRCspp_tradeoff_all, marker='o', color='blue', markersize = 6, linewidth = 3)
plt.plot('Tradeoff', 'Fnuc', data=CRCspp_tradeoff_all, marker='o', color='darkviolet', markersize = 6, linewidth = 3)
plt.plot('Tradeoff', 'Pana', data=CRCspp_tradeoff_all, marker='o', color='crimson', markersize = 6, linewidth = 3)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5)) #this moves legend to outside of the graph
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.xlabel("Tradeoff value", labelpad = 15)
plt.ylabel("Growth rate", labelpad = 15)
plt.title("MICOM- predicted individual and community growth \n rates for five species of CRC-associated microbes", size=BIGGER_SIZE)
#once things are ready, do this:
plt.savefig("/Users/adamo010/Documents/microbiome_FBA/CRC_indiv_plus_comm_GRs_01.29.20.png", dpi=300, bbox_inches="tight", pad_inches=0.1)


######################################################################################
######################################################################################
######################################################################################
#another thing I would like to do is see if I can compare the growth media of different
#communities. Do inputs/outputs differ? Goal: make a graph of the differences 

#try healthy vs CRC at community GRs 1.0, 0.8, and 0.6. 
#redo step 1 above

healthyspp_files = {"Fprau": "/Users/adamo010/Documents/microbiome_FBA/Fprau.json", 
              "Bbif": "/Users/adamo010/Documents/microbiome_FBA/Bbif.json",
              "Ccat": "/Users/adamo010/Documents/microbiome_FBA/Ccat.json",
              "Lacid": "/Users/adamo010/Documents/microbiome_FBA/Lacid.json",
              "Rinul": "/Users/adamo010/Documents/microbiome_FBA/Rinul.json"}









