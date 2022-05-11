#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 10:29:36 2020

@author: adamo010
"""
#update from 02.10.20: now have a contaminant-filtered OTU table: see 03.04.20_OTU_table_filtering_pipeline_clean.py
#update from 03.05.20: now have an otu table where XXX_level_sort columns don't have extra '; ' substrings when classification is incomplete

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

#But right now Spyder isn't working so will probably have to trash my computer anyway. 
#fixed. Seriously, I hate virtual environments

agora = pd.read_csv("AGORA_model_phylogeny_table.csv")
otus = pd.read_csv("03.10.20_Burns_2015_OTUS_contaminants_removed.csv") #this is the new contaminant-free table

#I think a new plan is in order

###############step 1: add columns to agora that can be matched to otus
#using the columns created in the contaminant filtering pipeline
adapted_agora = copy.deepcopy(agora)  #copy to keep everything clean

adapted_agora["phylum_level_sort"] = "k__" + adapted_agora['kingdom'] + ";  p__" + adapted_agora['phylum'] #add an extra column to get a k__ in front of kingdom
adapted_agora["class_level_sort"] = adapted_agora["phylum_level_sort"] + ";  c__" + adapted_agora['class']
adapted_agora["order_level_sort"] = adapted_agora["class_level_sort"] + ";  o__" + adapted_agora['order']
adapted_agora["family_level_sort"] = adapted_agora["order_level_sort"] + ";  f__" + adapted_agora['family']
adapted_agora["genus_level_sort"] = adapted_agora["family_level_sort"] + ";  g__" + adapted_agora['genus']
adapted_agora["species_level_sort"] = adapted_agora["genus_level_sort"] + ";  s__" + adapted_agora['species']
#add genus2 in there- don't know if that will help at all
adapted_agora["genus2_level_sort"] = adapted_agora["family_level_sort"] + ";  g__" + adapted_agora['genus2']
adapted_agora["species2_level_sort"] = adapted_agora["genus2_level_sort"] + ";  s__" + adapted_agora['species']

#check that this won't generate '; ;' issues like in otu table
for elem in adapted_agora['species_level_sort']:
    if elem[-1] == ";":
        print(elem)
#all good

###############step 2: create some dictionaries of model_id and XX_level_sort: SCRATCHED
#agora_species2_dict = adapted_agora.set_index('model_ID').to_dict()['species2_level_sort']
#agora_species_dict = adapted_agora.set_index('model_ID').to_dict()['species_level_sort']
#agora_genus2_dict = adapted_agora.set_index('model_ID').to_dict()['genus2_level_sort']
#agora_genus_dict = adapted_agora.set_index('model_ID').to_dict()['genus_level_sort']
#agora_family_dict = adapted_agora.set_index('model_ID').to_dict()['family_level_sort']
#each dictionary should have 818 key-value pairs- this is how many models are in the agora database. Redundancy at each taxonomic level is
#expected at this point (hence why we have to create lists of models)

###############step 3: create dictionaries for each taxonomic level
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

#####ALL RIGHT. What do we have now.
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
family_matched_models.to_csv("03.18.20_spp_genus_fam_level_models_mapped.csv")
family_matched_models.to_csv("03.18.20_spp2_genus2_fam2_level_models_mapped.csv")

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

V1_list = [otus_with_species_models, otus_with_genus_models, otus_with_family_models]
otus_with_models_V1 = pd.concat(V1_list)
otus_with_models_V1.to_csv("03.18.20_OTUs_with_models_V1.csv")

V2_list = [otus_with_species2_models, otus_with_genus2_models, otus_with_family2_models]
otus_with_models_V2 = pd.concat(V2_list)
otus_with_models_V2.to_csv("03.18.20_OTUs_with_models_V2.csv")


#without models is going to be more complicated- probably just pull the family and family2 lists
otus_without_family_models = family_matched_models[family_matched_models['family_matched_model_list'].isnull()]  
otus_without_family2_models = family2_matched_models[family2_matched_models['family_matched_model_list'].isnull()]   

#collect all other otus
otus_with_low_tax_ids = otus[otus['family'].isnull() & otus['genus'].isnull() & otus['species'].isnull()]   

wo_models = [otus_without_family_models, otus_with_low_tax_ids]
otus_without_models = pd.concat(wo_models)
otus_without_models.to_csv("03.18.20_OTUs_without_models_V1.csv")
wo_models2 = [otus_without_family2_models, otus_with_low_tax_ids]
otus_without_models2 = pd.concat(wo_models2)
otus_without_models2.to_csv("03.18.20_OTUs_without_models_V2.csv")











########################################old stuff from 02.10.20

#step 1: add a new column of empty lists to otus UPDATE DON'T ACTUALLY NEED THIS
#otus['model_list'] = np.empty((len(otus), 0)).tolist()
#for row in otus['model_list']:
    #print(type(row))

#step 2: add models to the model list YIKES

agora_genus_dict = agora.set_index('genus').to_dict()['model_ID']
#this creates a dictionary where the genus is the key and model_ID is the value
#NOPE this won't work, because every key has to be different and there are duplicate genuses.
agora_genus_dict2 = agora.set_index('model_ID').to_dict()['genus']
#THIS creates a dictionary where the model_ID is the value and genus is the key
#hopefully this will make the lookup easier.    
#I don't know if this is useful

#maybe I can pre-create genus lists of models. If genus is shared, append model_id to the genus list
unique_agora_genera = agora.genus.unique()

unique_agora_genera2 = pd.DataFrame(unique_agora_genera)    #import to pandas
unique_agora_genera2['model_list'] = np.empty((len(unique_agora_genera2), 0)).tolist()  #add an empty 'model list' column
unique_agora_genera2.reset_index()
unique_agora_genera2= unique_agora_genera2.rename(columns={0: "genus"})       #rename column

#what I want to do:
#for each row in unique_agora_genera2:
#scan through the genus column in agora
#if the value in the genus column in agora matches the value in the genus column of unique_agora_genera2,
#append the value of model_id in the agora dataframe to the list in the model_list column associated with that genus

models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]

for row in unique_agora_genera2['genus']:  
    genus_name_key = str(row)
    temp_df = agora[agora["genus"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    #now, I want to pull out the model names from these and append them to unique_agora_genera2       
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    models_dict[genus_name_key] = templist
    del templist    
    del genus_name_key

#EAT MY SHORTS PYTHON it finally worked. Precious models_dict
#now, do a better job of specifying different levels
genus_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_agora_genera2['genus']:  
    genus_name_key = str(row)
    temp_df = agora[agora["genus"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    #now, I want to pull out the model names from these and append them to unique_agora_genera2       
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    genus_models_dict[genus_name_key] = templist
    del templist    
    del genus_name_key
            
genus2_models_dict ={}         #this is where I'll store the key-value pairs for genus- [list of models]
for row in unique_agora_genera2['genus']:  
    genus_name_key = str(row)
    temp_df = agora[agora["genus2"].str.contains(str(row))]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    #now, I want to pull out the model names from these and append them to unique_agora_genera2       
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    del temp_df
    genus2_models_dict[genus_name_key] = templist
    del templist    
    del genus_name_key
                
#interesting- looks like genus vs genus2 (AGORA-exact or renamed to match new species names, respectively) don't matter that much
    
#now for the horrid part: matching to the otu data in otus

model_matched_otus = copy.deepcopy(otus)    

model_matched_otus['genus_matched_model_list'] = model_matched_otus['genus'].map(genus_models_dict)        
        
#actually, no problem because I'm great.
#let's take a look in excel because I'm actually not that great
model_matched_otus.to_csv("Burns_2015_OTUs_matched_models.csv")

############################################################################################
#OK, it's a new day and I have a new plan
#need to do this iteratively:
#match models to taxa with species-level data
#pull those taxa out of the OTUs list
#pull out taxa WITH species-level information that LACK a model- may make models for these
#match models to (remaning) taxa with genus-level data (lists of matching models required)
#pull those taxa out of the OTUs list
#pull out taxa WITH genus-level information that LACK a model- may make models for these
#match models to (remaning) taxa with family-level data (lists of matching models required)
#pull those taxa out of the OTUs list
#pull out taxa WITH family level information that LACK a model- may make models for these
#see where we're at as far as coverage

#part 0: remove all obvious contaminants
#see OTU_table_filtering_pipeline_development.py






#################### JUNK
for row in unique_agora_genera2['genus']:
    #print(row)
    if str(row) in agora['genus']:
        print(row)



for row in otus['genus']:         #for each row in the column taxon_ID in the dataframe otus:
   print(row)
   if row =  
    
    
    model_list_dict.update({str(row),[]})

#testing     
    temp_df = agora[agora["genus"].str.contains('Aeromonas')]        #temp_df should now contain only agora rows that match that row in unique_agora_genera2 
    #now, I want to pull out the model names from these and append them to unique_agora_genera2       
    templist = []
    for row in temp_df['model_ID']:
        templist.append(str(row))
    models_dict['Aeromonas'] = templist
    del temp_df
    del templist    
        


for index, row in filtered_otus.iterrows():   
    col_list = colnames
    counter2 = 0
    if col in colnames:
        print(col)
        for col in col_list:
            if row >= 0.001:
                counter2 = counter2+1
            fitered_otus["num_above_cutoff"] = counter2    
    
#LowAb_otus = preprocessed_otus[preprocessed_otus]        
#preprocessed_otus.drop(LowAb_otus, inplace=True)   

#for col in preprocessed_otus.columns:
    #counter = 0
    #if 'Sample_' in str(col):
       # for row in col:
          #  if row >=0.001:
           #     counter= counter+1
    i#f counter >=18:
#

#Ways NOT to remove blank taxonomic levels
otus_blanked.loc[len(otus_blanked['species'])>=4, 'species'] == ""

if len(otus_blanked['species']) <=4:
    otus_blanked.replace("s__ ", "")

otus_blanked.species = otus_blanked.species.replace({"s__": 0})
otus_blanked.species.replace("s__", "", inplace=True)

otus_blanked.loc[(otus_blanked.species == "s__"), "species"] = ''
otus_blanked['species'].mask(otus_blanked['species'] == "s__", "", inplace=True)



######## DON'T NEED THIS ANYMORE
#now, we have several dictionaries (unique_agora_species/species2/genus/genus2/family) and otus
#split otus into different dataframes identified to different taxonomic levels
#first, delete all cells with p__, c__, o__, f__, g__, s__ (recall that genus and spp have two columns)
#assign models as possible from AGORA, see if other databases have missing models.    

for elem in otus['species']:
    print(elem)
    print(len(elem))
#so, if the element length is less than or equal to 4, it's blank (essentially)

otus_blanked =copy.deepcopy(otus)

#another option NONE OF THESE WORK
#so I genuinely do not understand this, but just print all values less than four characters, save that
#value, and replace it. 
for elem in otus_blanked['species']:
    if len(elem) <= 4:
        #print(elem)
        dummy = elem
otus_blanked['species'].replace({dummy: ""}, inplace = True)        

def remove_blanks(otu_table, colname):
    for elem in otu_table[str(colname)]:
        if len(elem) <=4:
            dummy = elem
    otu_table[str(colname)].replace({dummy: ""}, inplace = True)
    del dummy
    return

remove_blanks(otus_blanked, "kingdom")
remove_blanks(otus_blanked, "phylum")
remove_blanks(otus_blanked, "class")
remove_blanks(otus_blanked, "order")
remove_blanks(otus_blanked, "family")
remove_blanks(otus_blanked, "genus")
remove_blanks(otus_blanked, "species")    

#okay, this worked really well. I'm wondering if it belongs earlier in the process (i.e. in OTU_table_filtering_pipeline_clean)

if pd.notna(otus['species']):
    print(row)
    
#otus_to_species_level = otus.groupby(pd.notna(otus['species']))





    
    