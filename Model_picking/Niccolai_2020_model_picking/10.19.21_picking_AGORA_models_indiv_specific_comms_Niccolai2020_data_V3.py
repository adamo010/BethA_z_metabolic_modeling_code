#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 11:13:02 2021

@author: adamo010
"""

#New in 10.19 version- a new contaminant table
#NEW in 10.15 version- a new contaminant table
#this is based on 07.13.21_Picking_AGORA_models_indiv_specific_comms_Hale2018_data.py. This version is for Niccolai2020 data
#this is based on 07.28.20_Picking_AGORA_models_indov_specific_comms.py. This version is for Hale 2018 data, instead of Burns 2015 data. 

#this is based on 06.16.20_Picking_AGORA_models_round_1.py. I've revised the project to focus on individual-specific microbiome communities,
#rather than an 'average' microbiome community for CRC vs healthy patients. 


#SOME HELPFUL NOTES:
#-we're using the MICOM_CRC_FBA_02.2020 virtual environment
#-the best home folder is /Users/adamo010/Documents/MICOM_CRC_FBA

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

#step 1: import data.
os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Niccolai_2020_model_picking")

agora = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/AGORA_model_phylogeny_table.csv")
otus = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Contaminant_removal/10.19.21_Niccolai2020_otu_table_abundance_filtered_0.1_any_sample_cutoff.csv") #this is the new contaminant-free table

###############step 1: add columns to agora that can be matched to otus
#using the columns created in the contaminant filtering pipeline
adapted_agora = copy.deepcopy(agora)  #copy to keep everything clean

#quick programming note: spacing changed here. earlier copies of this code have two spaces in the '; p' etc. scripts.
#one difference for this dataset is that there are NO species in the Niccolai OTU table. So we don't need "species level sort" stuff
adapted_agora["phylum_level_sort"] = adapted_agora['kingdom'] + "; " + adapted_agora['phylum'] #add an extra column to get a k__ in front of kingdom
adapted_agora["class_level_sort"] = adapted_agora["phylum_level_sort"] + "; " + adapted_agora['class']
adapted_agora["order_level_sort"] = adapted_agora["class_level_sort"] + "; " + adapted_agora['order']
adapted_agora["family_level_sort"] = adapted_agora["order_level_sort"] + "; " + adapted_agora['family']
adapted_agora["genus_level_sort"] = adapted_agora["family_level_sort"] + "; " + adapted_agora['genus']
#adapted_agora["species_level_sort"] = adapted_agora["genus_level_sort"] + "; " + adapted_agora['species']
#add genus2 in there- don't know if that will help at all
adapted_agora["genus2_level_sort"] = adapted_agora["family_level_sort"] + "; " + adapted_agora['genus2']
#adapted_agora["species2_level_sort"] = adapted_agora["genus2_level_sort"] + "; " + adapted_agora['species']

####intermediate step from 03.26.20: strip whitespace and see if this helps with model picking issue
adapted_agora['phylum_level_sort'] = adapted_agora['phylum_level_sort'].str.strip()
adapted_agora['class_level_sort'] = adapted_agora['class_level_sort'].str.strip()
adapted_agora['order_level_sort'] = adapted_agora['order_level_sort'].str.strip()
adapted_agora['family_level_sort'] = adapted_agora['family_level_sort'].str.strip()
adapted_agora['genus_level_sort'] = adapted_agora['genus_level_sort'].str.strip()
#adapted_agora['species_level_sort'] = adapted_agora['species_level_sort'].str.strip()
adapted_agora['genus2_level_sort'] = adapted_agora['genus2_level_sort'].str.strip()
#adapted_agora['species2_level_sort'] = adapted_agora['species2_level_sort'].str.strip()

#check that this won't generate '; ;' issues like in otu table
for elem in adapted_agora['genus_level_sort']:
    if elem[-1] == ";":
        print(elem)
del elem        
#all good
adapted_agora.to_csv("adapted_agora.csv")

###############step 2: create dictionaries for each taxonomic level
#e.g. key is a specific, say, genus; key is a list of models corresponding to that genus
#however, it is going to be most useful, I think, to make each key an XX_level_sort value
#that avoids the issue of different phyla containing genera with the same name
#again, not doing species level here. 

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
    del templist, genus_name_key, row
del unique_agora_genus

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
    del templist, genus2_name_key, row
del unique_agora_genus2

#family level
unique_agora_family = adapted_agora.family_level_sort.unique()
#something is wrong here: Bacillales Incertae Sedis XI is not being counted as unique
#never mind, Gemella had its order AND family changed
#adapted_agora.to_csv("adapted_agora.csv")
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
    del templist, family_name_key, row
del unique_agora_family

###############step 3: create new otu tables based on 'id to species level', 'id to genus level', 'id to family level' without replacement
#side note: have to do some extra steps because why wouldn't I.   
#NOTE that this code is a bit different from the 07.28.20 version because my OTU table looks different. Also, i'm lazy, so I'll just merge in
#some stuff from other dataframes

#import OTUs with taxonomy information
taxonomy_table = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Niccolai2020/Niccolai_2020_OTUs_with_matched_taxonomy.csv")
taxonomy_table =taxonomy_table.drop(columns={"Unnamed: 0"})

#NOTE: this table, for this dataset, DOES include read counts. We don't want these. (Relative abundances are already in the otus dataframe)
#keep first column and last 11 columns.
taxonomy_table.drop(taxonomy_table.columns[1:81], axis=1, inplace=True)

#nice. now we can merge
otus_plus_tax = pd.merge(left= otus, right= taxonomy_table, how="left", left_on="OTU_ID", right_on="OTU_ID")

#now, create new dataframes where otus are identified to the genus and family levels ONLY (i.e. get rid of low tax ids)
#otus_to_species_level = otus_plus_tax[otus_plus_tax['Species'].notnull()]   #hey, look at that
#now for the harder part
otus_to_genus_level = otus_plus_tax[otus_plus_tax['Genus'].notnull()]   #you may award me the genuis medal now
otus_to_family_level = otus_plus_tax[otus_plus_tax['Family'].notnull() & otus_plus_tax['Genus'].isnull()]   

#now, each of these dataframes needs to have a new column added where their xxx_level_sort column is matched to a value from xxx_models_dict
#fortunately, old me already figured this one out

#I thought we would only deepcopy the species2/ genus2, but because of a 'settingWithCopyWarning', I'll copy everything.
#swomething about otus_to_xxx_level above being slices of otus, rather than separate dataframes. Python is weird. 
#species_matched_models = copy.deepcopy(otus_to_species_level)
#species2_matched_models = copy.deepcopy(otus_to_species_level)
genus_matched_models = copy.deepcopy(otus_to_genus_level)
genus2_matched_models = copy.deepcopy(otus_to_genus_level)
family_matched_models = copy.deepcopy(otus_to_family_level)

#step 0: remove white space. Don't know if this is needed, but maybe?
#species_matched_models['species_level_sort'] = species_matched_models['species_level_sort'].str.strip()
#species_matched_models['genus_level_sort'] = species_matched_models['genus_level_sort'].str.strip()
#species_matched_models['family_level_sort'] = species_matched_models['family_level_sort'].str.strip()

#species2_matched_models['species_level_sort'] = species2_matched_models['species_level_sort'].str.strip()
#species2_matched_models['genus_level_sort'] = species2_matched_models['genus_level_sort'].str.strip()
#species2_matched_models['family_level_sort'] = species2_matched_models['family_level_sort'].str.strip()

#genus_matched_models['species_level_sort'] = genus_matched_models['species_level_sort'].str.strip()
genus_matched_models['genus_level_sort'] = genus_matched_models['genus_level_sort'].str.strip()
genus_matched_models['family_level_sort'] = genus_matched_models['family_level_sort'].str.strip()

#genus2_matched_models['species_level_sort'] = genus2_matched_models['species_level_sort'].str.strip()
genus2_matched_models['genus_level_sort'] = genus2_matched_models['genus_level_sort'].str.strip()
genus2_matched_models['family_level_sort'] = genus2_matched_models['family_level_sort'].str.strip()

#family_matched_models['species_level_sort'] = family_matched_models['species_level_sort'].str.strip()
family_matched_models['genus_level_sort'] = family_matched_models['genus_level_sort'].str.strip()
family_matched_models['family_level_sort'] = family_matched_models['family_level_sort'].str.strip()

####step 1: map species models to species-level list
#all right, here is where we're going to need to be careful with Niccolai2020 data. There's some weird shit going on at the genus level.
#that means there are multiple species and genus ids. How do we handle this?
#fortunately, Hale2018 taught us how to EXPLODE
#EXPLODE!!!
#first, convert Genus column to list (instead of strings delimited by /)
def convert_to_list(x):
	return x.split("/")
genus_matched_models["Genus"] = genus_matched_models["Genus"].apply(convert_to_list)
genus2_matched_models["Genus"] = genus2_matched_models["Genus"].apply(convert_to_list)
#gross, for some reason it only works if I write a function to split the string, and then apply it to the column in question

#now explode!
genus_matched_models= genus_matched_models.explode('Genus')
genus2_matched_models= genus2_matched_models.explode('Genus')

#spectacular. Now, back to matching models.
#whew, this might take some figuring, given that we don't have species. That's okay. 
#step 0: new species_level_sort column for this. Do I... need this for genus level sort? Maybe not. B/c we're not matching "genus-species"
#genus_matched_models["clean_genus_level_sort"] = genus_matched_models["genus_level_sort"] + "; " + species_matched_models['Species']
#genus2_matched_models["clean_genus_level_sort"] = genus2_matched_models["genus_level_sort"] + "; " + species2_matched_models['Species']

####step 1: map genus models to genus-level list
genus_matched_models['genus_matched_model_list'] = genus_matched_models['genus_level_sort'].map(genus_models_dict)        
genus2_matched_models['genus2_matched_model_list'] = genus2_matched_models['genus_level_sort'].map(genus2_models_dict)        

####step 2:  pull out otus which are identified to the genus level but which lack a genus-level model
no_genus_models = genus_matched_models[genus_matched_models['genus_matched_model_list'].isnull()]
no_genus2_models = genus2_matched_models[genus2_matched_models['genus2_matched_model_list'].isnull()]   

####step 3: append modelless speces-level otus to genus dataframe
#genus_matched_models = genus_matched_models.append(no_species_models, ignore_index = True)
#genus2_matched_models = genus2_matched_models.append(no_species2_models, ignore_index = True)

####step 4: map genus models to genus-level (plus modelless species-level) list
#genus_matched_models['genus_matched_model_list'] = genus_matched_models['genus_level_sort'].map(genus_models_dict)        
#genus2_matched_models['genus2_matched_model_list'] = genus2_matched_models['genus_level_sort'].map(genus2_models_dict)      

####step 5: pull out otus which are identified to the genus level but which lack a genus-level model
#no_genus_models = genus_matched_models[genus_matched_models['genus_matched_model_list'].isnull()]   
#no_genus2_models = genus2_matched_models[genus2_matched_models['genus2_matched_model_list'].isnull()]   

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
#otus_with_species_models = species_matched_models[species_matched_models['species_matched_model_list'].notnull()]   
#otus_with_species2_models = species2_matched_models[species2_matched_models['species2_matched_model_list'].notnull()]   
otus_with_genus_models = genus_matched_models[genus_matched_models['genus_matched_model_list'].notnull()]   
otus_with_genus2_models = genus2_matched_models[genus2_matched_models['genus2_matched_model_list'].notnull()]   
otus_with_family_models = family_matched_models[family_matched_models['family_matched_model_list'].notnull()]   
otus_with_family2_models = family2_matched_models[family2_matched_models['family_matched_model_list'].notnull()]   

#without models is going to be more complicated- probably just pull the family and family2 lists
otus_without_family_models = family_matched_models[family_matched_models['family_matched_model_list'].isnull()]  
otus_without_family2_models = family2_matched_models[family2_matched_models['family_matched_model_list'].isnull()]   
#collect all other otus
otus_with_low_tax_ids = otus_plus_tax[otus_plus_tax['Family'].isnull() & otus_plus_tax['Genus'].isnull()]   

V1_list = [otus_with_genus_models, otus_with_family_models]
otus_with_models_V1 = pd.concat(V1_list)
V2_list = [otus_with_genus2_models, otus_with_family2_models]
otus_with_models_V2 = pd.concat(V2_list)

#####make some files to look at.
otus_with_models_V1.to_csv("10.19.21_Niccolai2020_OTUs_with_models_V1_indiv_specific_comms_V3.csv")
otus_with_models_V2.to_csv("10.19.21_Niccolai2020_OTUs_with_models_V2_indiv_specific_comms_V3.csv")
otus_with_low_tax_ids.to_csv("10.19.21_Niccolai2020_OTUs_with_low_tax_ids_V3.csv")
otus_without_family_models.to_csv("10.19.21_Niccolai2020_OTUs_without_models_V1_indiv_specific_comms_V3.csv") 
otus_without_family2_models.to_csv("10.19.21_Niccolai2020_OTUs_without_models_V2_indiv_specific_comms_V3.csv") 

#let's take a look and see what we get!


