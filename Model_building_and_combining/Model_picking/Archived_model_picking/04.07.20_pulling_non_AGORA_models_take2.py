#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 16:29:32 2020

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


#04.07.20 UPDATE: have gone in by hand and extracted all model phylogenies that I can. By googling them, essentially.
#note that I have extracted and new species names- phylogenies have been updated. Will probably just end up using the old species
#names, since my dataset is also old. But I have the code to do multiple species, so why the hell not

#note that the next step after this is to build models, so enjoy the 'I've figured this shit out already' while you can.

#We'll be pulling from Building_microbiome_communities_03.25.20.py

mambo = pd.read_csv("04.09.20_MAMBO_model_spp_for_id.csv")       #found a bunch of taxonomy issues, so I fixed it today.
otus = pd.read_csv("03.27.20_OTUs_without_models_V1.csv") #at least V1 and V2 are identical

###############step 1: add columns to mambo that can be matched to otus
#using the columns created in the contaminant filtering pipeline
adapted_mambo = copy.deepcopy(mambo)  #copy to keep everything clean

adapted_mambo["phylum_level_sort"] = "k__" + adapted_mambo['Kingdom'] + ";  p__" + adapted_mambo['Phylum'] #add an extra column to get a k__ in front of kingdom
adapted_mambo["class_level_sort"] = adapted_mambo["phylum_level_sort"] + ";  c__" + adapted_mambo['Class']
adapted_mambo["order_level_sort"] = adapted_mambo["class_level_sort"] + ";  o__" + adapted_mambo['Order']
adapted_mambo["family_level_sort"] = adapted_mambo["order_level_sort"] + ";  f__" + adapted_mambo['Family']
adapted_mambo["genus_level_sort"] = adapted_mambo["family_level_sort"] + ";  g__" + adapted_mambo['Genus_extracted']
adapted_mambo["genus2_level_sort"] = adapted_mambo["family_level_sort"] + ";  g__" + adapted_mambo['Genus_new']
adapted_mambo["species_level_sort"] = adapted_mambo["genus_level_sort"] + ";  s__" + adapted_mambo['Species_extracted']
adapted_mambo["species2_level_sort"] = adapted_mambo["genus2_level_sort"] + ";  s__" + adapted_mambo['Species_extracted']

#important note: genus/species is the genus/species extracted from the model itself
#genus2/species2 is the latest taxonomic assignment according to NCBI

####intermediate step: strip whitespace and see if this helps with model picking issue
adapted_mambo['phylum_level_sort'] = adapted_mambo['phylum_level_sort'].str.strip()
adapted_mambo['class_level_sort'] = adapted_mambo['class_level_sort'].str.strip()
adapted_mambo['order_level_sort'] = adapted_mambo['order_level_sort'].str.strip()
adapted_mambo['family_level_sort'] = adapted_mambo['family_level_sort'].str.strip()
adapted_mambo['genus_level_sort'] = adapted_mambo['genus_level_sort'].str.strip()
adapted_mambo['species_level_sort'] = adapted_mambo['species_level_sort'].str.strip()
adapted_mambo['genus2_level_sort'] = adapted_mambo['genus2_level_sort'].str.strip()
adapted_mambo['species2_level_sort'] = adapted_mambo['species2_level_sort'].str.strip()

#need to add a new thing that deletes NaNs and leaves blank strings. Not really sure why this suddenly came up, but here we are.
adapted_mambo = adapted_mambo.replace(np.nan, "", regex=True)

###############step 2: create dictionaries for each taxonomic level
#e.g. key is a specific, say, genus; key is a list of models corresponding to that genus
#however, it is going to be most useful, I think, to make each key an XX_level_sort value
#that avoids the issue of different phyla containing genera with the same name

unique_mambo_species = adapted_mambo.species_level_sort.unique()
#then, do other stuff
unique_mambo_species = pd.DataFrame(unique_mambo_species)    #import to pandas
unique_mambo_species['model_list'] = np.empty((len(unique_mambo_species), 0)).tolist()  #add an empty 'model list' column
unique_mambo_species.reset_index()
unique_mambo_species= unique_mambo_species.rename(columns={0: "species_level_sort"})       #rename column

#then, create a dictionary where each key is a unique species and each value is a list of models that match that species
species_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_mambo_species['species_level_sort']:  
    species_name_key = str(row)
    temp_df = adapted_mambo[adapted_mambo["species_level_sort"].str.contains(str(row))]        #temp_df should now contain only mambo rows that match that row in unique_mambo_genera2 
    #now, I want to pull out the model names from these and append them to unique_mambo_genera2       
    templist = []
    for row in temp_df['model_file_name']:
        templist.append(str(row))
    del temp_df
    species_models_dict[species_name_key] = templist
    del templist    
    del species_name_key

#could I put this into a function? Yes. Would it take longer than copy-pasting? Also yes. Am I going to be doing this repeatedly enough
#that a function would be useful? No
unique_mambo_species2 = adapted_mambo.species2_level_sort.unique()
#then, do other stuff
unique_mambo_species2 = pd.DataFrame(unique_mambo_species2)    #import to pandas
unique_mambo_species2['model_list'] = np.empty((len(unique_mambo_species2), 0)).tolist()  #add an empty 'model list' column
unique_mambo_species2.reset_index()
unique_mambo_species2= unique_mambo_species2.rename(columns={0: "species2_level_sort"})       #rename column
#then, create a dictionary where each key is a unique species and each value is a list of models that match that species
species2_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]

for row in unique_mambo_species2['species2_level_sort']:  
    species2_name_key = str(row)
    temp_df = adapted_mambo[adapted_mambo["species2_level_sort"].str.contains(str(row))]        #temp_df should now contain only mambo rows that match that row in unique_mambo_genera2 
    #now, I want to pull out the model names from these and append them to unique_mambo_genera2       
    templist = []
    for row in temp_df['model_file_name']:
        templist.append(str(row))
    del temp_df
    species2_models_dict[species2_name_key] = templist
    del templist    
    del species2_name_key


##############MIGHT WORK LATER- working on stuff above right now

#genus level
unique_mambo_genus = adapted_mambo.genus_level_sort.unique()
#then, do other stuff
unique_mambo_genus = pd.DataFrame(unique_mambo_genus)    #import to pandas
unique_mambo_genus['model_list'] = np.empty((len(unique_mambo_genus), 0)).tolist()  #add an empty 'model list' column
unique_mambo_genus.reset_index()
unique_mambo_genus= unique_mambo_genus.rename(columns={0: "genus_level_sort"})       #rename column

#then, create a dictionary where each key is a unique genus and each value is a list of models that match that genus
genus_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_mambo_genus['genus_level_sort']:  
    genus_name_key = str(row)
    temp_df = adapted_mambo[adapted_mambo["genus_level_sort"].str.contains(str(row))]        #temp_df should now contain only mambo rows that match that row in unique_mambo_genera2 
    #now, I want to pull out the model names from these and append them to unique_mambo_genera2       
    templist = []
    for row in temp_df['model_file_name']:
        templist.append(str(row))
    del temp_df
    genus_models_dict[genus_name_key] = templist
    del templist    
    del genus_name_key

#genus2 level
unique_mambo_genus2 = adapted_mambo.genus2_level_sort.unique()
unique_mambo_genus2 = pd.DataFrame(unique_mambo_genus2)    #import to pandas
unique_mambo_genus2['model_list'] = np.empty((len(unique_mambo_genus2), 0)).tolist()  #add an empty 'model list' column
unique_mambo_genus2.reset_index()
unique_mambo_genus2= unique_mambo_genus2.rename(columns={0: "genus2_level_sort"})       #rename column
genus2_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_mambo_genus2['genus2_level_sort']:  
    genus2_name_key = str(row)
    temp_df = adapted_mambo[adapted_mambo["genus2_level_sort"].str.contains(str(row))]        #temp_df should now contain only mambo rows that match that row in unique_mambo_genera2 
    #now, I want to pull out the model names from these and append them to unique_mambo_genera2       
    templist = []
    for row in temp_df['model_file_name']:
        templist.append(str(row))
    del temp_df
    genus2_models_dict[genus2_name_key] = templist
    del templist    
    del genus2_name_key    

#family level
unique_mambo_family = adapted_mambo.family_level_sort.unique()
unique_mambo_family = pd.DataFrame(unique_mambo_family)    #import to pandas
unique_mambo_family['model_list'] = np.empty((len(unique_mambo_family), 0)).tolist()  #add an empty 'model list' column
unique_mambo_family.reset_index()
unique_mambo_family= unique_mambo_family.rename(columns={0: "family_level_sort"})       #rename column
family_models_dict ={}         #this is where I'll store the key-value pairs for family- [list of models]
for row in unique_mambo_family['family_level_sort']:  
    family_name_key = str(row)
    temp_df = adapted_mambo[adapted_mambo["family_level_sort"].str.contains(str(row))]        #temp_df should now contain only mambo rows that match that row in unique_mambo_genera2 
    #now, I want to pull out the model names from these and append them to unique_mambo_genera2       
    templist = []
    for row in temp_df['model_file_name']:
        templist.append(str(row))
    del temp_df
    family_models_dict[family_name_key] = templist
    del templist    
    del family_name_key

#peachy


#next step: create new otu tables based on 'id to species level', 'id to genus level', 'id to family level' without replacement

otus_to_species_level = otus[otus['species'].notna()]   #hey, look at that

#now for the harder part
otus_to_genus_level = otus[otus['genus'].notna() & otus['species'].isnull()]   #you may award me the genuis medal now
otus_to_family_level = otus[otus['family'].notna() & otus['genus'].isnull() & otus['species'].isnull()]   
#11/10

#now, each of these dataframes needs to have a new column added where their xxx_level_sort column is matched to a value from xxx_models_dict
#fortunately, old me already figured this one out

#I thought we would only deepcopy the species2/ genus2, but because of a 'settingWithCopyWarning', I'll copy everything.
#swomething about otus_to_xxx_level above being slices of otus, rather than separate dataframes. Python is weird. 
species_matched_models = copy.deepcopy(otus_to_species_level)
species2_matched_models = copy.deepcopy(otus_to_species_level)
genus_matched_models = copy.deepcopy(otus_to_genus_level)
genus2_matched_models = copy.deepcopy(otus_to_genus_level)
family_matched_models = copy.deepcopy(otus_to_family_level)

#step 0: remove white space
species_matched_models['species_level_sort'] = species_matched_models['species_level_sort'].str.strip()
species_matched_models['genus_level_sort'] = species_matched_models['genus_level_sort'].str.strip()
species_matched_models['family_level_sort'] = species_matched_models['family_level_sort'].str.strip()

species2_matched_models['species_level_sort'] = species2_matched_models['species_level_sort'].str.strip()
species2_matched_models['genus_level_sort'] = species2_matched_models['genus_level_sort'].str.strip()
species2_matched_models['family_level_sort'] = species2_matched_models['family_level_sort'].str.strip()

genus_matched_models['species_level_sort'] = genus_matched_models['species_level_sort'].str.strip()
genus_matched_models['genus_level_sort'] = genus_matched_models['genus_level_sort'].str.strip()
genus_matched_models['family_level_sort'] = genus_matched_models['family_level_sort'].str.strip()

genus2_matched_models['species_level_sort'] = genus2_matched_models['species_level_sort'].str.strip()
genus2_matched_models['genus_level_sort'] = genus2_matched_models['genus_level_sort'].str.strip()
genus2_matched_models['family_level_sort'] = genus2_matched_models['family_level_sort'].str.strip()

family_matched_models['species_level_sort'] = family_matched_models['species_level_sort'].str.strip()
family_matched_models['genus_level_sort'] = family_matched_models['genus_level_sort'].str.strip()
family_matched_models['family_level_sort'] = family_matched_models['family_level_sort'].str.strip()

####step 1: map species models to species-level list
species_matched_models['species_matched_model_list'] = species_matched_models['species_level_sort'].map(species_models_dict)        
species2_matched_models['species2_matched_model_list'] = species2_matched_models['species_level_sort'].map(species2_models_dict)        

####step 2: pull out otus which are identified to the species level but which lack a species-level model
no_species_models = species_matched_models[species_matched_models['species_matched_model_list'].isnull()]
no_species2_models = species2_matched_models[species2_matched_models['species2_matched_model_list'].isnull()]   

####step 3: append modelless speces-level otus to genus dataframe
genus_matched_models = genus_matched_models.append(no_species_models, ignore_index = True)
genus2_matched_models = genus2_matched_models.append(no_species2_models, ignore_index = True)

####step 4: map genus models to genus-level (plus modelless species-level) list
genus_matched_models['genus_matched_model_list'] = genus_matched_models['genus_level_sort'].map(genus_models_dict)        
genus2_matched_models['genus2_matched_model_list'] = genus2_matched_models['genus_level_sort'].map(genus2_models_dict)      

####step 5: pull out otus which are identified to the genus level but which lack a genus-level model
no_genus_models = genus_matched_models[genus_matched_models['genus_matched_model_list'].isnull()]   
no_genus2_models = genus2_matched_models[genus2_matched_models['genus2_matched_model_list'].isnull()]   

####step 6: append modelless genus-level otus to family dataframe
family_matched_models = family_matched_models.append(no_genus_models, ignore_index = True)
family2_matched_models = family_matched_models.append(no_genus2_models, ignore_index = True)

####step 7: prepare family level model lists (need 2, for genus/species and genus2/species2)
family2_matched_models = copy.deepcopy(family_matched_models)  

####step 8: map family models to family-level list (plus species-level and genus-level OTUs without models)
family_matched_models['family_matched_model_list'] = family_matched_models['family_level_sort'].map(family_models_dict)     
family2_matched_models['family_matched_model_list'] = family_matched_models['family_level_sort'].map(family_models_dict)      
#save these as CSVs: will want to look at later
family_matched_models.to_csv("03.26.20_spp_genus_fam_level_models_mapped.csv")
family_matched_models.to_csv("03.26.20_spp2_genus2_fam2_level_models_mapped.csv")

####step 9: pull out otus which are identified to the family level but which lack a family-level model
no_family_models = family_matched_models[family_matched_models['family_matched_model_list'].isnull()]   
no_family2_models = family2_matched_models[family2_matched_models['family_matched_model_list'].isnull()]   

#great, now we have some ots (at species/genus/family level) with OTUs assigned
#need to compile two lists- one of otus with modesl assigned, and one of otus without

#if species_matched_model_list column is not blank, append to new dataframe
otus_with_species_models = species_matched_models[species_matched_models['species_matched_model_list'].notnull()]   
otus_with_species2_models = species2_matched_models[species2_matched_models['species2_matched_model_list'].notnull()]   
otus_with_genus_models = genus_matched_models[genus_matched_models['genus_matched_model_list'].notnull()]   
otus_with_genus2_models = genus2_matched_models[genus2_matched_models['genus2_matched_model_list'].notnull()]   
otus_with_family_models = family_matched_models[family_matched_models['family_matched_model_list'].notnull()]   
otus_with_family2_models = family2_matched_models[family2_matched_models['family_matched_model_list'].notnull()]   

#to export to csvs: these CSVs are JUST mambo-matched models. Will append to agora later
V1_list = [otus_with_species_models, otus_with_genus_models, otus_with_family_models]
otus_with_models_V1 = pd.concat(V1_list)
otus_with_models_V1.to_csv("04_08.20_OTUs_with_mambo_models_V1.csv")

V2_list = [otus_with_species2_models, otus_with_genus2_models, otus_with_family2_models]
otus_with_models_V2 = pd.concat(V2_list)
otus_with_models_V2.to_csv("04_08.20_OTUs_with_mambo_models_V2.csv")

#without models is going to be more complicated- probably just pull the family and family2 lists
otus_without_family_models = family_matched_models[family_matched_models['family_matched_model_list'].isnull()]  
otus_without_family2_models = family2_matched_models[family2_matched_models['family_matched_model_list'].isnull()]   

#collect all other otus
otus_with_low_tax_ids = otus[otus['family'].isnull() & otus['genus'].isnull() & otus['species'].isnull()]   

wo_models = [otus_without_family_models, otus_with_low_tax_ids]
otus_without_models = pd.concat(wo_models)
otus_without_models.to_csv("04_08.20_OTUs_without_mambo_models_V1.csv")
wo_models2 = [otus_without_family2_models, otus_with_low_tax_ids]
otus_without_models2 = pd.concat(wo_models2)
otus_without_models2.to_csv("04_08.20_OTUs_without_mambo_models_V2.csv")

#cool. how many extra models did we get?    

