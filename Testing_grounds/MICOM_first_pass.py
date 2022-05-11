#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 16:21:06 2020

@author: adamo010
"""

#step 0: in terminal: pip install micom
#step 0.1: in terminal: install gurobi
#conda config --add channels http://conda.anaconda.org/gurobi
#conda install gurobi
#grbgetkey ed7e7426-3d57-11ea-a2f6-0a7c4f30bdbe
#(this is the academic license I signed up for on 01.22.20)
#note that I had to import it into the virtual environment to get it to run

#unfuckingbelievable. Shouldn't have to do this b/c gurobi should be part of Anaconda, but here we are.
#go to /Library/gurobi900/mac64 in the terminal (not in the virtual environment)
#type python setup.py install
#type import gurobipy from inside the reactivated virtual environment

import micom
import gurobipy
import scipy
import numpy as np

#Just gonna run through the tutorial on this one
#########step 1: building communities###########

from micom.data import test_taxonomy
input_spp = test_taxonomy()
input_spp
type(input_spp)
#test_taxonomy is a community generated from a test dataframe. The dataframe requires
#at least two columns: ID (which specifies the taxon ID) and file (which specifies where the model is)
#test_taxonomy contains five E.coli models with six features: id, genus, species, reactions,
#metabolites, and file. Reactions and metabolites are numberical values; the rest are string

from micom import Community #this function creates a community from a pandas series of models
fivespp = Community(input_spp)   #the community named fivespp will be made from the models in taxonomy by the Community function
print("Built a community with a total of {} reactions.".format(len(fivespp.reactions)))
#this prints a statement that identifies the number of reactions in com

#now, a community is built; however, the taxa in this community are still stored
fivespp.taxonomy
#this also contains a column called abundance; by default, each taxon is present at an identical abundance
#abundances can be specified in the original taxonomy if need be.

#we can also import models from AGORA:
from micom.data import agora
agoraspp = agora.copy()
agoraspp.file= "models/" + agoraspp.file
agoraspp.head()
#tax is now a pandas dataframe containing all the AGORA models

#it is computationally intensive to construct large community models; therefore,
#once they are generated, they should be saved in a serialized format (whatever that is)

%time fivespp = Community(input_spp)
#takes 1.34s to create a five species community library
%time fivespp.to_pickle("fivespp.pickle")
#save community as a pickle file: community.pickle: takes 14.4ms
from micom import load_pickle
%time fivespp = load_pickle("fivespp.pickle")
#now, when the pickle community file is loaded, takes 72.7ms (much shorter than assembling de novo)

#######step 2: community growth rates and fluxes########
from micom import Community, data

input_spp2 = data.test_taxonomy()
fivespp2 = Community(input_spp2, solver='gurobi')
#note that the tutorial said to use gurobi for the solver, but I didn't have it installed.
#wtf, why would it suddenly work now?

print(fivespp2.objective.expression)
#by default, the objective for a community model is the community growth rate

fivespp2.optimize()
#prints out the community model solution
#tutorial claims that it prints out a table of indiv spp, but I didn't see it

sol = fivespp2.optimize()
sol.members
#this DOES print out a table of community members, with their respective abundances,
#growth rates, Rxs, and metabolites

#MICOM does not do fluxes by default (too slow), but it CAN if you insist.
sol = fivespp2.optimize(fluxes=True)
sol.fluxes

#MICOM will also do parsimonious FBA if you like
sol = fivespp2.optimize(fluxes=True, pfba= True)
sol.fluxes

############step 3: cooperative tradeoff#########
#need gurobi for this so have to figure out how to load it in.
#this is where community growth rate is lowered to boost indiv growth rate

CTsol = fivespp2.cooperative_tradeoff(fraction=1.0)
#yeah, this doesnt work with glpk. but how do I set gurobi as the default?
#everything in this stupid tutorial fails the first time, then I mess around with it, then it works the first way I tried it.
CTsol

#set some growth rate limits
CTsol_min1 = fivespp2.cooperative_tradeoff(min_growth=0.1)  # single value
CTsol_min1
CTsolmin2 = fivespp2.cooperative_tradeoff(min_growth=[0.1, 0.2, 0.3, 0.4, 0.5])
CTsolmin2

#test the impact of the tradeoff parameter on your solution
#make sure numpy is engaged
CTsol_tradeoff = fivespp2.cooperative_tradeoff(fraction=np.arange(0.1, 1.01, 0.1))
CTsol_tradeoff
#THIS IS THE COOLEST FEATURE

#"the solution can then be inspected by the usual pandas means" that everyone knows duh
rates = CTsol_tradeoff.solution.apply(lambda x: x.members.growth_rate)
rates
#what this does is goes through the solution column in healthyspp_CT_tradeoff and iterates through each column.
#for each value of solution, apply the lambda function (x represents the solution): get the members.growth rate
#another way to do this:
for row in rates['solution']:
    print(row)
    #print(row.members)      #prints each community members' abundance, growth rate, Rxs, and metabolites
    print(row.members.growth_rate)  #prints each community members' growth rate only 
    #print(row.members.abundance)
    #print(row.members.reactions)
    #print(row.members.metabolites)

#NEAT


############step 4: growth media##############
from micom import Community, data

input_spp3 = data.test_taxonomy()
fivespp3 = Community(input_spp3, solver='gurobi')

#analyse the minimal growth medium by enforcing a particular community growth rate (here, 0.8)
from micom.media import minimal_medium
fivesppMM = minimal_medium(fivespp3, 0.8)
fivesppMM
#this solution gives the smallest total import flux

#to minimize the number of import fluxes used:
minimal_medium(fivespp3, 0.8, minimize_components=True)
###this one isn't working the way it shoud on the tutorial...

#anyway. Cooperative tradeoff results can also be used as constraints in the MM calculation
fivesppMM_CT = fivespp3.cooperative_tradeoff()
growthrates = fivesppMM_CT.members.growth_rate.drop("medium")  # extracellular medium has no growth rate
fivesppMM_CT_med = minimal_medium(fivespp3, fivesppMM_CT.growth_rate, min_growth=growthrates)
fivesppMM_CT_med

#getting the minimal medium when lowering the required growth rates (here to 95% of optimum)
med95 = minimal_medium(fivespp3, 0.95*fivesppMM_CT.growth_rate, min_growth=0.95*growthrates)
med95

#obtaining export fluxes: negative values are exports, positive values are imports
medout = minimal_medium(fivespp3, fivesppMM_CT.growth_rate, min_growth=growthrates, exports=True)
medout

#applying growth media: use a dictionary structure that maps exchange Rx to their upper import flux bound (e.g. the output of minimal_medium)
inputmed = minimal_medium(fivespp3, 0.8, min_growth=0.8)   
#create a dictionary called inputmed, which contains the output of minimal_medium with the community fivespp3, 
#a community growth rate of 0.8, and a minimum growth of 0.8
fivespp3.medium = inputmed
#set the medium of the community fivespp3 to be inputmed
fivespp3.optimize()
#the output of fivespp3 shows a growth rate (0.8) that matches the GR defined when we calculated inputmed

###########step 5: taxa knockouts#########
#goal: identify cooperative and competitive interactions among spp by removing each one
#if a KO of one spp increases the GR of another, they are competing; if it decreases the GR, they are cooperating
#for communities of 100-500 taxa, full-factorial KOs can take 1-5h.

from micom import Community, data
input_spp4 = data.test_taxonomy()
fivespp4 = Community(input_spp4, solver="gurobi")

ko = fivespp4.knockout_species(fraction=1.0)
ko 
#by default, knockout_species returns a change in growth rate (before-after KO GRs)
#rows indicate the taxon knocked out, and columns indicate the change in GR

ko2 = fivespp4.knockout_species(fraction=1.0, diag=False)
ko2 
#to suppress the diagonal entries (which are zeroes anyway), add the diag=False parameter

ko_raw = fivespp4.knockout_species(fraction=1.0, method="raw")
ko_raw 
#this gives the straight-up new growth rates after the KO

ko_relative = fivespp4.knockout_species(fraction=1.0, method="relative change")
ko_relative 
#this gives the relative change ((new-old)/old)

#########step 6: messing with the models########
##6a: what impact does changing the abundance have on exchanges?


##6b: what impact does changing the diet (exchange bounds) have on exchanges
from micom import Community, data, elasticity
from micom import elasticity
from micom.elasticity import elasticities  # this is a dimensionless, normalized measure of how much a given parameter affects a flux value
#note that the tutorial states exchange_elasticities shoud be imported, but this isn't in the code at all

input_spp5 = data.test_taxonomy()
fivespp5 = Community(input_spp5, solver="gurobi")

eps = elasticities(fivespp5, fraction=1.0)
eps.head()
#this doesn't look like it's supposed to, but I'm not sure how to follow the manual when exchange_elasticities isn't a defined function












