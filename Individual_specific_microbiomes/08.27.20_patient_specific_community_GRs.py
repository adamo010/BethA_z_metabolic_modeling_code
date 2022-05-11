#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 11:25:18 2020

@author: adamo010
"""

#Have made community pickle files. Time to run some GROWTH RATES
#use the directory /Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes

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


########write some pilot code and see how this goes.
#Step 1: import pickle files (3)
from micom import load_pickle

pickle_file_list = []

#for file in os.listdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Pickle_sample_folder"):
    #src= file
    #if fnmatch.fnmatch(file, "*_community.pickle"):
        #pickle_file_list.append(str(file))

#pickle_file_list now contains all the pickle files
sample_name_list = []
for elem in pickle_file_list:
    elem2 = elem.rstrip("_community.pickle")
    sample_name_list.append(elem2)        

pickle_file_dict = dict(zip(sample_name_list, pickle_file_list))

#import pickle files: started at 9:55am
model_list = []
for key, value in pickle_file_dict.items():
    model = load_pickle(value)
    model_list.append(model)
    
#############start here after pickles import

model_list_dict = dict(zip(sample_name_list, model_list))    

fba_output_list = []

for key, value in model_list_dict.items():
    sol = value.optimize()  
    fba_output_list.append(sol)
    
fba_output_dict = dict(zip(sample_name_list, fba_output_list))   

#############################REAL THING########################

#Step 1: import pickle files
from micom import load_pickle

#get all .pickle files
pickle_file_list = []

Sample_11_comm = load_pickle("Sample_11_community.pickle")
Sample_12_comm = load_pickle("Sample_12_community.pickle")
Sample_13_comm = load_pickle("Sample_13_community.pickle")

pickle_file_list.append(Sample_11_comm)
pickle_file_list.append(Sample_12_comm)
pickle_file_list.append(Sample_13_comm)

pickle_file_names = ["Sample_11_comm", "Sample_12_comm", "Sample_13_comm"]

pickle_file_solutions =[]
for elem in pickle_file_list:
    sol = elem.optimize()
    pickle_file_solutions.append(sol)
    
fba_output_dict = dict(zip(pickle_file_names, pickle_file_solutions))   

for key, value in fba_output_dict.items():
    print(value.members)

#now I would like to print this....    
#was going to do this in a full on "for key, value in..." thing, but I really think a function is the way to go.
#I did it in a function in microbiome_fba_on_diets.py and it worked, and doing it differently here is NOT working, so...
    
#the format will be as follows:
#for key, value in fba_output_dict.items():
    #model_modifications(value)    
#value in these cases is the fba output     
    
#def model_modifications(fba_output_value, fba_output_key):
    #output = open("_prefba_output_summary.csv", 'w')
   # output.write(str(fba_output_value.members))
   # output.close()
   # with open("_fba_output_summary.csv", "w") as csvFile:         #create a new Diet_flux_summary file with columns aligned properly etc
        #writer = csv.writer(csvFile)
        #with open("_prefba_output_summary.csv") as input:
       #     fluxreader = csv.reader(input) 
            #for row in fluxreader:
               # row = row[0].split()
               # firstel = row[0]
               # if len(row) < 5 and firstel == "nan":
                  #  row.insert(0, "nan")
                  #  print(row)
               # writer.writerow(row)                    #create a new Diet_flux_summary file with columns aligned properly etc
   # os.remove("_prefba_output_summary.csv")      
   # for file in os.listdir():                       
     #   src=file
      #  if fnmatch.fnmatch(file, "_fba_output_summary.csv"):
      #      dst = str(fba_output_key)+file
       #     os.rename(src,dst)
    #return        


#for key, value in fba_output_dict.items():
   # model_modifications(value)    

#############
#still getting the stupid IndexError: list index out of rante on the row = row[0].split() command. I think it's because there's a blank row
#FFS.... okay, let's try it with ONE FUCKING SAMPLE FILE  
Sample_11_solution = Sample_11_comm.optimize()   
 
output = open("_prefba_output_summary.csv", 'w')
results = Sample_11_solution.members #this is a pandas dataframe
pd.set_option('display.max_rows', None)
output.write(str(results))
output.close()

with open("_fba_output_summary.csv", "w") as csvFile:         #create a new Diet_flux_summary file with columns aligned properly etc
    writer = csv.writer(csvFile)
    with open("_prefba_output_summary.csv") as input:
        fluxreader = csv.reader(input) 
        print(fluxreader)
        fluxreader2 = filter(None, fluxreader) #get rid of empty rows
        for row in fluxreader2:
            row = row[0].split()
            firstel = row[0]
            if len(row) < 5:
                row = ['compartments', 'abundance', 'growth_rate', 'reactions', 'metabolites']
            writer.writerow(row)                    #create a new Diet_flux_summary file with columns aligned properly etc
df = pd.read_csv('_fba_output_summary.csv', skiprows=1)
df.to_csv('_fba_output_summary.csv', index=False)
os.remove("_prefba_output_summary.csv")      
for file in os.listdir():                       
    src=file
    if fnmatch.fnmatch(file, "_fba_output_summary.csv"):
        dst = "Sample_11"+file
        os.rename(src,dst)
    
#well, that's a start. But now I'm getting the weird file cutoff thing where the middle part of the output is just "..."
        #Fixed! Thanks pd.set_option('display.max_rows', None)
     
######### NOW, let us get back to doing this in parallel. Want to write a function to run this for every solution in fba_output_dict
        
def model_modifications(fba_output_value, fba_output_key):
    output = open("_prefba_output_summary.csv", 'w')
    results = fba_output_value.members #this is a pandas dataframe
    pd.set_option('display.max_rows', None)
    output.write(str(results))
    output.close()
    del results
    del output
    with open("_fba_output_summary.csv", "w") as csvFile:         #create a new Diet_flux_summary file with columns aligned properly etc
        writer = csv.writer(csvFile)
        with open("_prefba_output_summary.csv") as input:
            fluxreader = csv.reader(input) 
            fluxreader2 = filter(None, fluxreader) #get rid of empty rows
            for row in fluxreader2:
                row = row[0].split()
                firstel = row[0]
                if len(row) < 5:
                    row = ['compartments', 'abundance', 'growth_rate', 'reactions', 'metabolites']
                writer.writerow(row)                    #create a new Diet_flux_summary file with columns aligned properly etc
    df = pd.read_csv('_fba_output_summary.csv', skiprows=1)
    df.to_csv('_fba_output_summary.csv', index=False)
    os.remove("_prefba_output_summary.csv")      
    for file in os.listdir():                       
        src=file
        if fnmatch.fnmatch(file, "_fba_output_summary.csv"):
            dst = str(fba_output_key)+file
            os.rename(src,dst)
    del src
    del dst
    del file
    del row
    del fluxreader
    del fluxreader2
    del firstel
    return        

for key, value in fba_output_dict.items():
    model_modifications(value, key)    
 
#well, that... kind of worked. need to adjust the names on those columns, but it's actually not too bad. 
    #Donezo!
    
#######    

fba_output_file_list = []
for file in os.listdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes"):
    src= file
    if fnmatch.fnmatch(file, "*comm_fba_output_summary.csv"):
        fba_output_file_list.append(str(file))
            
fbaed_sample_name_list = []
for elem in fba_output_file_list:
    elem2 = elem.rstrip("_fba_output_summary.csv")
    fbaed_sample_name_list.append(elem2)        

fba_output_file_dict = dict(zip(fbaed_sample_name_list, fba_output_file_list))

fba_output_pd_series = []
for key, value in fba_output_file_dict.items():
    pdseries = pd.read_csv(value)
    fba_output_pd_series.append(pdseries)
    
fba_output_pd_dict = dict(zip(fbaed_sample_name_list, fba_output_pd_series))

#I think I'm really screwing myself with these file naming conventions    

#okay, where is my dictionary of otus (keys) and taxonomies (values)?
preprocessed_otus = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal/Burns_absolute_OTU_table.csv") #note that this is in /Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal
preprocessed_otus.rename(columns={"#OTU_ID":"OTU_ID"}, inplace = True)
#only keep first and last column
otu_taxonomy_lookup_table = copy.deepcopy(preprocessed_otus)
otu_taxonomy_lookup_table = otu_taxonomy_lookup_table[['OTU_ID','taxonomy']]
#need to convert this to string for downstream stuff
otu_taxonomy_lookup_table["OTU_ID"] = otu_taxonomy_lookup_table["OTU_ID"].astype(str)

#now, I will use otu_taxonomy_lookup_table to match taxonomies to the otus in the fba_output_pd_series dataframes
#first, i need to make sure the fba_output_pd_series have the correct column headers

for key, value in fba_output_pd_dict.items():
    value.columns=["OTU_ID", "abundance", "growth_rate", "reactions", "metabolites"]
    print(value.columns)
       
#df1.merge(df2, left_on='lkey', right_on='rkey')

#for key, value in fba_output_pd_dict.items():
    otu_table = otu_taxonomy_lookup_table.OTU_ID.astype(str)
    output = pd.merge(left = value, right= otu_table, left_on ="OTU_ID", right_on ="OTU_ID")
    print(output)
    output.to_csv("_fba_results_plus_taxonomies.csv", index=False)
    for file in os.listdir():                       
        src=file
        if fnmatch.fnmatch(file, "_fba_results_plus_taxonomies.csv"):
            dst = str(key)+file
            os.rename(src,dst)

for key, value in fba_output_pd_dict.items():
    output = pd.merge(value, otu_taxonomy_lookup_table, on="OTU_ID", how="left")
    #print(output)
    output.to_csv("_fba_results_plus_taxonomies.csv", index=False)
    for file in os.listdir():                       
        src=file
        if fnmatch.fnmatch(file, "_fba_results_plus_taxonomies.csv"):
            dst = str(key)+file
            os.rename(src,dst)

#WORKS #SCIENCE FRIDAY
       
#########################    JUNK     #########################    JUNK     #########################    JUNK

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

###########

#that actually worked pretty nicely.
#next, I would like to tack on a column to each of these files where the taxonomy is included (right now it's just OTU ids)
#have to import the csv, match up a dictionary of otu_id/taxonomy, then re-save the csv files

fba_output_file_list = []
for file in os.listdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes"):
    src= file
    if fnmatch.fnmatch(file, "*_fba_output_summary.csv"):
        fba_output_file_list.append(str(file))
            
fbaed_sample_name_list = []
for elem in fba_output_file_list:
    elem2 = elem.rstrip("_fba_output_summary.csv")
    fbaed_sample_name_list.append(elem2)        

fba_output_file_dict = dict(zip(fbaed_sample_name_list, fba_output_file_list))

fba_output_pd_series = []
for key, value in fba_output_file_dict.items():
    pdseries = pd.read_csv(value)
    fba_output_pd_series.append(pdseries)
    
fba_output_pd_dict = dict(zip(fbaed_sample_name_list, fba_output_pd_series))

#I think I'm really screwing myself with these file naming conventions    

#okay, where is my dictionary of otus (keys) and taxonomies (values)?
preprocessed_otus = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal/Burns_absolute_OTU_table.csv") #note that this is in /Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal
preprocessed_otus.rename(columns={"#OTU_ID":"OTU_ID"}, inplace = True)
#only keep first and last column
otu_taxonomy_lookup_table = preprocessed_otus[['OTU_ID','taxonomy']]

#now, I will use otu_taxonomy_lookup_table to match taxonomies to the otus in the fba_output_pd_series dataframes
#first, i need to make sure the fba_output_pd_series have the correct colum headers

for key, value in fba_output_pd_dict.items():
    value.columns=["OTU_ID", "abundance", "growth_rate", "reactions", "metabolites"]
    for col in value.columns:
        print(col[0])
    
    
#df1.merge(df2, left_on='lkey', right_on='rkey')

for key, value in fba_output_pd_series.items():
    merged_pd = value.merge(left = fba_output_pd_series, right= otu_taxonomy_lookup_table, how="left", left_on ="" )


    

##########################################################JUNK

for file in os.listdir():
    src= file
    if fnmatch.fnmatch(file, "*_community.pickle"):
        pickle_file_list.append(str(file))
        
#pickle_file_list now contains all the pickle files
sample_name_list = []
for elem in pickle_file_list:
    elem2 = elem.rstrip("_community.pickle")
    sample_name_list.append(elem2)        

#zip into a dictionary
pickle_file_dict = dict(zip(sample_name_list, pickle_file_list))

#import pickle files
for key, value in pickle_file_dict.items():
    key = load_pickle(value)
#took two days.

print(elem2)

#sol = com.optimize()    

    
import sys
original_stdout = sys.stdout
    
with open("_prefba_output_summary.csv", 'w') as f:
    sys.stdout = f
    print(str(Sample_11_solution.members))
    sys.stdout = original_stdout