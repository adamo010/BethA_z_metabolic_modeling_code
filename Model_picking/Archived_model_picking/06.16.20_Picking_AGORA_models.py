#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 10:29:36 2020

@author: adamo010
"""

#DO NOT USE ME!!!!!! USe 06.16.20_picking_AGORA_models_round2_fixed_taxonomies.
#update from 04.27.20: changed contaminants list, so ended up with slightly different OTU input table. 

#update from 04.27.20: NO NEW MODELS PICKED so no need to change downstream pipeline.
#just want to get a csv file of all the otus that didn't ahve sufficiently deep taxonomic id to have models picked/
#combined/built for them. Found in step 3.

#update fro 04.21.20: freshly edited contaminant list, and square brackets removed. 
#update from 04.15.20: I don't honestly know. Just trying to keep the dates congruent.
#okay, I have it- remove 03.26.20 updates. Keep it simple, stupid. with only 180 OTUs, can fix them by hand if needed.
#update from 04.13.20: using an abundance and contaminant-filtered OTU table.
#update from 03.25.20: used to be called Building_microbiome_communities_03.25.20.py. Note also that this adds in 
#Building_microbiome_communities_take2_03.26.20.py
#update from 03.10.20: stupid fucking code is missing a bunch of models for some stupid ass reason that I don't understand.
#update from 03.05.20: now have an otu table where XXX_level_sort columns don't have extra '; ' substrings when classification is incomplete
#update from 02.10.20: now have a contaminant-filtered OTU table: see 03.04.20_OTU_table_filtering_pipeline_clean.py

#anyway, to install gurobi:
#step 0: in terminal: pip install micom
#step 0.1: in terminal: install gurobi
#conda config --add channels http://conda.anaconda.org/gurobi
#conda install gurobi
#grbgetkey ed7e7426-3d57-11ea-a2f6-0a7c4f30bdbe
#(this is the academic license I signed up for on 01.22.20)
#note that I had to import it into the virtual environment to get it to run

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

#now that we have MICOM running and some real data, need to figure out how to build communities.
#I think the pipeline will be as follows:
    #import pandas series of phylogenetic breakdown from Burns OTU table and AGORA model list
    #add a new column in Burns OTU table where each cell is a(n empty for now) list
    #for each (start with genus, then go to family) in the Burns OTU table, find all models with a matching genus in
    #the AGORA model list, and append those models to the new list column in the Burns OTU table

agora = pd.read_csv("AGORA_model_phylogeny_table.csv")
otus = pd.read_csv("06.16.20_otu_table_abundance_filtered_0.1.csv") #this is the new contaminant-free table

#remove square brackets from all taxon ids. In excel. Because it's a clusterfuck in python.

###############step 1: add columns to agora that can be matched to otus
#using the columns created in the contaminant filtering pipeline
adapted_agora = copy.deepcopy(agora)  #copy to keep everything clean

#quick programming note: spacing changed here. earlier copies of this code have two spaces in the '; p' etc. scripts.
adapted_agora["phylum_level_sort"] = "k__" + adapted_agora['kingdom'] + "; p__" + adapted_agora['phylum'] #add an extra column to get a k__ in front of kingdom
adapted_agora["class_level_sort"] = adapted_agora["phylum_level_sort"] + "; c__" + adapted_agora['class']
adapted_agora["order_level_sort"] = adapted_agora["class_level_sort"] + "; o__" + adapted_agora['order']
adapted_agora["family_level_sort"] = adapted_agora["order_level_sort"] + "; f__" + adapted_agora['family']
adapted_agora["genus_level_sort"] = adapted_agora["family_level_sort"] + "; g__" + adapted_agora['genus']
adapted_agora["species_level_sort"] = adapted_agora["genus_level_sort"] + "; s__" + adapted_agora['species']
#add genus2 in there- don't know if that will help at all
adapted_agora["genus2_level_sort"] = adapted_agora["family_level_sort"] + "; g__" + adapted_agora['genus2']
adapted_agora["species2_level_sort"] = adapted_agora["genus2_level_sort"] + "; s__" + adapted_agora['species']

####intermediate step from 03.26.20: strip whitespace and see if this helps with model picking issue
adapted_agora['phylum_level_sort'] = adapted_agora['phylum_level_sort'].str.strip()
adapted_agora['class_level_sort'] = adapted_agora['class_level_sort'].str.strip()
adapted_agora['order_level_sort'] = adapted_agora['order_level_sort'].str.strip()
adapted_agora['family_level_sort'] = adapted_agora['family_level_sort'].str.strip()
adapted_agora['genus_level_sort'] = adapted_agora['genus_level_sort'].str.strip()
adapted_agora['species_level_sort'] = adapted_agora['species_level_sort'].str.strip()
adapted_agora['genus2_level_sort'] = adapted_agora['genus2_level_sort'].str.strip()
adapted_agora['species2_level_sort'] = adapted_agora['species2_level_sort'].str.strip()

#check that this won't generate '; ;' issues like in otu table
for elem in adapted_agora['species_level_sort']:
    if elem[-1] == ";":
        print(elem)
del elem        
#all good

###############step 2: create dictionaries for each taxonomic level
#e.g. key is a specific, say, genus; key is a list of models corresponding to that genus
#however, it is going to be most useful, I think, to make each key an XX_level_sort value
#that avoids the issue of different phyla containing genera with the same name

#first, create a unique taxa list: this one is for species
unique_agora_species = adapted_agora.species_level_sort.unique()
#then, do other stuff
unique_agora_species = pd.DataFrame(unique_agora_species)    #import to pandas
unique_agora_species['model_list'] = np.empty((len(unique_agora_species), 0)).tolist()  #add an empty 'model list' column
unique_agora_species.reset_index()
unique_agora_species= unique_agora_species.rename(columns={0: "species_level_sort"})       #rename column

#then, create a dictionary where each key is a unique species and each value is a list of models that match that species
species_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_agora_species['species_level_sort']:  
    species_name_key = str(row)
    temp_df = adapted_agora[adapted_agora["species_level_sort"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    #now, I want to pull out the model names from these and append them to unique_agora_genera2       
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    species_models_dict[species_name_key] = templist
    del templist    
    del species_name_key
    del row

#could I put this into a function? Yes. Would it take longer than copy-pasting? Also yes. Am I going to be doing this repeatedly enough
#that a function would be useful? No
unique_agora_species2 = adapted_agora.species2_level_sort.unique()
#then, do other stuff
unique_agora_species2 = pd.DataFrame(unique_agora_species2)    #import to pandas
unique_agora_species2['model_list'] = np.empty((len(unique_agora_species2), 0)).tolist()  #add an empty 'model list' column
unique_agora_species2.reset_index()
unique_agora_species2= unique_agora_species2.rename(columns={0: "species2_level_sort"})       #rename column
#then, create a dictionary where each key is a unique species and each value is a list of models that match that species
species2_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_agora_species2['species2_level_sort']:  
    species2_name_key = str(row)
    temp_df = adapted_agora[adapted_agora["species2_level_sort"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    #now, I want to pull out the model names from these and append them to unique_agora_genera2       
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    species2_models_dict[species2_name_key] = templist
    del templist    
    del species2_name_key
    del row

#genus level
unique_agora_genus = adapted_agora.genus_level_sort.unique()
#then, do other stuff
unique_agora_genus = pd.DataFrame(unique_agora_genus)    #import to pandas
unique_agora_genus['model_list'] = np.empty((len(unique_agora_genus), 0)).tolist()  #add an empty 'model list' column
unique_agora_genus.reset_index()
unique_agora_genus= unique_agora_genus.rename(columns={0: "genus_level_sort"})       #rename column

#then, create a dictionary where each key is a unique genus and each value is a list of models that match that genus
genus_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_agora_genus['genus_level_sort']:  
    genus_name_key = str(row)
    temp_df = adapted_agora[adapted_agora["genus_level_sort"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    #now, I want to pull out the model names from these and append them to unique_agora_genera2       
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    genus_models_dict[genus_name_key] = templist
    del templist    
    del genus_name_key
    del row

#genus2 level
unique_agora_genus2 = adapted_agora.genus2_level_sort.unique()
unique_agora_genus2 = pd.DataFrame(unique_agora_genus2)    #import to pandas
unique_agora_genus2['model_list'] = np.empty((len(unique_agora_genus2), 0)).tolist()  #add an empty 'model list' column
unique_agora_genus2.reset_index()
unique_agora_genus2= unique_agora_genus2.rename(columns={0: "genus2_level_sort"})       #rename column
genus2_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_agora_genus2['genus2_level_sort']:  
    genus2_name_key = str(row)
    temp_df = adapted_agora[adapted_agora["genus2_level_sort"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    #now, I want to pull out the model names from these and append them to unique_agora_genera2       
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    genus2_models_dict[genus2_name_key] = templist
    del templist    
    del genus2_name_key
    del row    

#family level
unique_agora_family = adapted_agora.family_level_sort.unique()
unique_agora_family = pd.DataFrame(unique_agora_family)    #import to pandas
unique_agora_family['model_list'] = np.empty((len(unique_agora_family), 0)).tolist()  #add an empty 'model list' column
unique_agora_family.reset_index()
unique_agora_family= unique_agora_family.rename(columns={0: "family_level_sort"})       #rename column
family_models_dict ={}         #this is where I'll store the key-value pairs for family- [list of models]
for row in unique_agora_family['family_level_sort']:  
    family_name_key = str(row)
    temp_df = adapted_agora[adapted_agora["family_level_sort"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    #now, I want to pull out the model names from these and append them to unique_agora_genera2       
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    family_models_dict[family_name_key] = templist
    del templist    
    del family_name_key
    del row

###############step 3: create new otu tables based on 'id to species level', 'id to genus level', 'id to family level' without replacement

otus_to_species_level = otus[otus['species'].notna()]   #hey, look at that

#now for the harder part
otus_to_genus_level = otus[otus['genus'].notna() & otus['species'].isnull()]   #you may award me the genuis medal now
otus_to_family_level = otus[otus['family'].notna() & otus['genus'].isnull() & otus['species'].isnull()]   

#############NEW as of 06.16.20: create a csv of all otus not identified to at least the family level
otus_lacking_fine_taxids = otus[otus['family'].isnull() & otus['genus'].isnull() & otus['species'].isnull()]   
otus_lacking_fine_taxids.to_csv("06.16.20_OTUs_lacking_fine_taxonomic_ids.csv")
#########end

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

####step 6: prepare family level model lists (need 2, for genus/species and genus2/species2)
family2_matched_models = copy.deepcopy(family_matched_models)  

####step 7: append modelless genus-level otus to family dataframe
family_matched_models = family_matched_models.append(no_genus_models, ignore_index = True)
family2_matched_models = family2_matched_models.append(no_genus2_models, ignore_index = True)

####step 8: map family models to family-level list (plus species-level and genus-level OTUs without models)
family_matched_models['family_matched_model_list'] = family_matched_models['family_level_sort'].map(family_models_dict)     
family2_matched_models['family_matched_model_list'] = family2_matched_models['family_level_sort'].map(family_models_dict)      
#save these as CSVs: will want to look at later- DO NOT DO THIS YET have take2 to go through
#family_matched_models.to_csv("03.26.20_spp_genus_fam_level_models_mapped.csv")
#family_matched_models.to_csv("03.26.20_spp2_genus2_fam2_level_models_mapped.csv")

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

#without models is going to be more complicated- probably just pull the family and family2 lists
otus_without_family_models = family_matched_models[family_matched_models['family_matched_model_list'].isnull()]  
otus_without_family2_models = family2_matched_models[family2_matched_models['family_matched_model_list'].isnull()]   
#collect all other otus
otus_with_low_tax_ids = otus[otus['family'].isnull() & otus['genus'].isnull() & otus['species'].isnull()]   

V1_list = [otus_with_species_models, otus_with_genus_models, otus_with_family_models]
otus_with_models_V1 = pd.concat(V1_list)
V2_list = [otus_with_species2_models, otus_with_genus2_models, otus_with_family2_models]
otus_with_models_V2 = pd.concat(V2_list)

#####make some files to look at.
otus_with_models_V1.to_csv("06.16.20_OTUs_with_models_V1.csv")
otus_with_models_V2.to_csv("06.16.20_OTUs_with_models_V2.csv")
otus_with_low_tax_ids.to_csv("06.16.20_OTUs_with_low_tax_ids.csv")
otus_without_family_models.to_csv("06.16.20_OTUs_without_models_V1.csv")
otus_without_family2_models.to_csv("06.16.20_OTUs_without_models_V2.csv")


