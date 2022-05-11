#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 17:30:02 2021

@author: adamo010
"""

#NEW in 10.19 version- a new contaminant table
#NEW in 10.15 version- a new contaminant table
#this is based on 07.15.21_Picking_MAMBO_models_indiv_specific_comms_Hale2018_data.py. This version is for Niccolai2020 data
#this is based on 07.28.20_Picking_MAMBO_models_indov_specific_comms.py. This version is for Hale 2018 data, instead of Burns 2015 data. 

#this is based on 06.16.20_Picking_MAMBO_models.py. I've revised the project to focus on individual-specific microbiome communities,
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

mambo = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/07.23.20_MAMBO_model_taxonomy_assigned_plus_Treponema.csv")
otus = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Model_picking/Niccolai_2020_model_picking/10.19.21_Niccolai2020_OTUs_without_models_V2_indiv_specific_comms_V3.csv") #this is the new contaminant-free table
otus.drop(columns={"Unnamed: 0", "taxonomy_orig_y", "keep_or_toss"}, inplace=True) #clean up this mess a bit

###############step 1: add columns to mambo that can be matched to otus
#using the columns created in the contaminant filtering pipeline
adapted_mambo = copy.deepcopy(mambo)  #copy to keep everything clean

#quick programming note: spacing changed here. earlier copies of this code have two spaces in the '; p' etc. scripts.
#one difference for this dataset is that there are NO species in the Niccolai OTU table. So we don't need "species level sort" stuff
adapted_mambo["phylum_level_sort"] = adapted_mambo['Kingdom'] + "; " + adapted_mambo['Phylum'] #add an extra column to get a k__ in front of kingdom
adapted_mambo["class_level_sort"] = adapted_mambo["phylum_level_sort"] + "; " + adapted_mambo['Class']
adapted_mambo["order_level_sort"] = adapted_mambo["class_level_sort"] + "; " + adapted_mambo['Order']
adapted_mambo["family_level_sort"] = adapted_mambo["order_level_sort"] + "; " + adapted_mambo['Family']
adapted_mambo["genus_level_sort"] = adapted_mambo["family_level_sort"] + "; " + adapted_mambo['Genus_extracted']
adapted_mambo["genus2_level_sort"] = adapted_mambo["family_level_sort"] + "; " + adapted_mambo['Genus_new']
#adapted_mambo["species_level_sort"] = adapted_mambo["genus_level_sort"] + "; " + adapted_mambo['Species_extracted']
#adapted_mambo["species2_level_sort"] = adapted_mambo["genus2_level_sort"] + "; " + adapted_mambo['Species_extracted']

####intermediate step: strip whitespace and see if this helps with model picking issue
adapted_mambo['phylum_level_sort'] = adapted_mambo['phylum_level_sort'].str.strip()
adapted_mambo['class_level_sort'] = adapted_mambo['class_level_sort'].str.strip()
adapted_mambo['order_level_sort'] = adapted_mambo['order_level_sort'].str.strip()
adapted_mambo['family_level_sort'] = adapted_mambo['family_level_sort'].str.strip()
adapted_mambo['genus_level_sort'] = adapted_mambo['genus_level_sort'].str.strip()
#adapted_mambo['species_level_sort'] = adapted_mambo['species_level_sort'].str.strip()
adapted_mambo['genus2_level_sort'] = adapted_mambo['genus2_level_sort'].str.strip()
#adapted_mambo['species2_level_sort'] = adapted_mambo['species2_level_sort'].str.strip()

#need to add a new thing that deletes NaNs and leaves blank strings. Not really sure why this suddenly came up, but here we are.
adapted_mambo = adapted_mambo.replace(np.nan, "", regex=True)

###############step 2: create dictionaries for each taxonomic level
#e.g. key is a specific, say, genus; key is a list of models corresponding to that genus
#however, it is going to be most useful, I think, to make each key an XX_level_sort value
#that avoids the issue of different phyla containing genera with the same name
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
    del templist, genus_name_key, row
del unique_mambo_genus    

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
    del templist, genus2_name_key, row
del unique_mambo_genus2 

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
    del templist, family_name_key, row
del unique_mambo_family

###############step 3: create new otu tables based on 'id to species level', 'id to genus level', 'id to family level' without replacement
#side note: have to do some extra steps because why wouldn't I.   
#NOTE that this code is a bit different from the 07.28.20 version because my OTU table looks different. Also, i'm lazy, so I'll just merge in
#some stuff from other dataframes

###############DOUBLE NOTE that, as long as you're using the output from picking AGORA models, taxonomy has already been merged in. 
#import OTUs with taxonomy information
#taxonomy_table = pd.read_csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Published_data/CRC_OTU_tables_Hale2018/Hale_2018_OTUs_with_matched_taxonomy.csv")
#taxonomy_table =taxonomy_table.drop(columns={"Unnamed: 0"})
#otus_plus_tax = pd.merge(left= otus, right= taxonomy_table, how="left", left_on="OTU_ID", right_on="OTU_ID")
################
#instead, copy OTUs to create otus_plus_tax
otus_plus_tax = copy.deepcopy(otus)

#NEW for 07.15.21 version MAMBO models: clean up this dataframe to get rid of extra columns.
otus_plus_tax = otus_plus_tax.drop(columns={"genus2_matched_model_list", "family_matched_model_list"})

#now, create new dataframes where otus are identified to the species, genus, and family levels ONLY (i.e. get rid of low tax ids)
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
otus_with_models_V1.to_csv("10.19.21_Niccolai2020_OTUs_with_MAMBO_models_V1_indiv_specific_comms_V3.csv")
otus_with_models_V2.to_csv("10.19.21_Niccolai2020_OTUs_with_MAMBO_models_V2_indiv_specific_comms_V3.csv")
otus_with_low_tax_ids.to_csv("10.19.21_Niccolai2020_OTUs_with_MAMBO_low_tax_ids_V3.csv")
otus_without_family_models.to_csv("10.19.21_Niccolai2020_OTUs_without_MAMBO_models_V1_indiv_specific_comms_V3.csv") 
otus_without_family2_models.to_csv("10.19.21_Niccolai2020_OTUs_without_MAMBO_models_V2_indiv_specific_comms_V3.csv") 

#let's take a look and see what we get!


