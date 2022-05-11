#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 10:06:57 2020

@author: adamo010
"""

#Okay, now I have a bunch of new CarveMe models and a bunch of smushed together models. Time to see if they actually run. 


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
import ast

#going to use cobra!
#getting a weird error when I try to import micom: AttributeError: module 'cobra.io.sbml' has no attribute 'LOGGER'
#maybe it'll be okay? I don't actually need MICOM right now. But this might be a concern later

from cobra.medium import minimal_medium #not sure if I'll need this but just in case

#goal is to basically run all the models in each folder and print out their results to a dictionary where I can check the results
#first, need to get a list of all the models in a given directory
from os import listdir
from os.path import isfile, join

#get all CarveMe models into a list called carved_models
onlyfiles = [f for f in listdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/0_genome_seqs_for_model_building") if isfile(join("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining/0_genome_seqs_for_model_building", f))]

carved_models = []
for filename in onlyfiles:
    if filename.endswith(".xml"):
        carved_models.append(filename)
del filename
del onlyfiles

#get all merged models into a list called merged_models
onlyfiles2 = [f for f in listdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining") if isfile(join("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_building_and_combining", f))]

merged_models = []
for filename in onlyfiles2:
    if filename.endswith("model.xml"):
        merged_models.append(filename)
del filename
del onlyfiles2

#cool. Now, I want to run FBA on all the models from carved_models and merged_models
#borrow some of the model-running code from microbiome_fba_cleancopy.py

model_test_results = {}

def metabolic_model(taxon_model):
    fba_model= cobra.io.read_sbml_model(taxon_model)
    fba_model_sol = fba_model.optimize()
    sol_string = str(fba_model_sol)
    model_test_results.update( {taxon_model : sol_string} )
    return

for model in merged_models:
    metabolic_model(model)
    

carved_model_test_results= {}

def metabolic_model2(taxon_model):
    fba_model= cobra.io.read_sbml_model(taxon_model)
    fba_model_sol = fba_model.optimize()
    sol_string = str(fba_model_sol)
    carved_model_test_results.update( {taxon_model : sol_string} )
    return

for model in carved_models:
    metabolic_model2(model)
    
#save models and take a look at the results
merged_model_results = open("05.26.20_merged_model_test_results.csv", "w")
writer= csv.writer(merged_model_results)
for key, value in model_test_results.items():
    writer.writerow([key, value])
    
merged_model_results.close()    
    
carved_model_results = open("05.26.20_carved_model_test_results.csv", "w")
writer= csv.writer(carved_model_results)
for key, value in carved_model_test_results.items():
    writer.writerow([key, value])
    
carved_model_results.close()   










def metabolic_model(taxon_model, modelname):
    #taxon_model is an .xml model
    #####part 1
    fba_model= cobra.io.read_sbml_model(taxon_model)      
    ########new addition: nested function for making maximal medium
    fba_model_sol = fba_model.slim_optimize() # run FBA and save the solution as model_sol
    print(fba_model.summary())
    #global fba_model_MM
    fba_model_MM = minimal_medium(fba_model, fba_model_sol)      #save the above solution (for the minimal medium needed for model to grow at its max growth rate) as the pandas series model_MM
    #global fba_model_MM_dict
    fba_model_MM_dict = fba_model_MM.to_dict()        #creates a dictionary of model's minimal Rxs and their fluxes
    ########
    diet = {}
    with open(diet_fluxes, "r") as file:
        next(file)
        for line in file:                                   
            linestuff = line.strip().split()                    
            diet[linestuff[0]] = float(linestuff[1])       #this creates a dictionary called 'diet' containing the metabolite (key) and its amount (value) from the supplied diet file
    import copy  
    diet_input = copy.deepcopy(diet)   #avoid editing the original file
    diet_input
    ######part 2
    fba_model_exchange_ids = [exchange.id for exchange in fba_model.exchanges]    #list comprehension: gets all the export Rxs (exchange ids) from the Faec_praus model
    #####part 3
    #at this point, we have diet (a list of input metabolites) and exchange IDs from the Faec_praus model. Need to reconcile them
    diet_input_matched = {}
    diet_input_matched = {k.replace("[e]", "(e)"): v for k, v in diet_input.items()}
    diet_input_matched
    #this fixes the problem that the diet list contains both [e] and (e), whereas the model exchange ids are only (e)
    diet_input_slimmed = {} 
    for key, value in diet_input_matched.items():                
        if key in fba_model_exchange_ids:          
            diet_input_slimmed[key] = value   
    diet_input_slimmed 
    #now we have a diet input for the model where the metabolites match what is present in the model
    missing_food = {}
    for key, value in fba_model_MM_dict.items():
        if key in diet_input_slimmed:
            print("ok")
        else:
            missing_food[key] = value
    missing_food  
    #this creates a dictionary of the metabolites and their amounts that ARE in the Faec_praus MM list, but are NOT in the diet input      
    diet_input_rounded = dict(diet_input_slimmed, **missing_food)  #nice
    #this crams the two lists together to make something that the model can actually run
    ##NEW EDIT: need to match amounts in diet_input_rounded to amounts in a minimal medium for 1g biomass
    fba_minimal = cobra.medium.minimal_medium(fba_model, 1, minimize_components = True) #thanks JC
    fba_minimal_dict = fba_minimal.to_dict()  
    #here, x is a dictionary = Fp_minimal_dict
    #here, y is a dictionary = diet_input_rounded
    x = fba_minimal_dict
    y = diet_input_rounded
    diet_input_min_growth = {key: max(x[key], value) if key in x.keys() else y[key] for key, value in y.items()}  
    #Jeremy saved my butt again: diet_input_min_growth contains all the metabolites and their amounts from
    #diet_input_rounded UNLESS the amount is lower than the amount in Fp_minimal; if this is the case, 
    #diet_input_min_growth contains the amount from Fp_minimal
    #I gotta learn list comprehension, seriously.
    fba_model.medium = diet_input_min_growth #no key error if everything worked
    #this is a good place to stop and print model_input_rounded. Should be a long list where the last 22ish items (after EX_2dmmq8) are constant between models
    fba_model.medium
    output = open("prediet_flux_summary.csv", "w")      #create a pre-diet file
    output.write(str(fba_model.summary()))             #write the model summmary to the prediet file
    output.close()                                      #close the prediet file
    with open("Diet_flux_summary.csv", "w") as csvFile:         #create a new Diet_flux_summary file with columns aligned properly etc
            writer = csv.writer(csvFile)
            with open("prediet_flux_summary.csv") as input:
                fluxreader = csv.reader(input) 
                for row in fluxreader:
                    row = row[0].split()
                    firstel = row[0]
                    if len(row) < 5 and firstel == "nan":
                        row.insert(0, "nan")
                        print(row)
                    writer.writerow(row)                    #create a new Diet_flux_summary file with columns aligned properly etc
    os.remove("prediet_flux_summary.csv")      
    fba_model_sol_diet = fba_model.optimize() # run FBA and save the solution as Faec_praus_sol_vegan
    fba_model_sol_diet.fluxes.to_csv("Diet_fluxes.csv", header=True)
    fba_model_sol_diet.shadow_prices.to_csv("Diet_shadow_prices.csv", header=True)
    fba_model_sol_diet.status
    for file in os.listdir():                       #rename files by their appropriate diet; couldn't figure out a way to do this in the above code
        src=file
        if fnmatch.fnmatch(file, "Diet_*.csv"):
            dst = str(sppname)+str(dietname)+file
            os.rename(src,dst)
            #shutil.move(source+f, dest)
    return
