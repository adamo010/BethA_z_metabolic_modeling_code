#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:52:11 2020

@author: adamo010
"""

#Okay, now we have three different rounds of model choosing that need to be concatenated:
#AGORA round1, then MAMBO, then AGORA round 2 with taxonomies fixed

#this code should also be useful when we generate other models and append them to different taxa.

#here we go!

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


agora1 = pd.read_csv("04.27.20_OTUs_with_models_V1.csv")
mambo = pd.read_csv("04.27.20_OTUs_with_mambo_models_V1.csv")
agora2 = pd.read_csv("05.01.20_OTUs_with_models_V1_taxonomy_fixed.csv")

otu_list = [agora1, mambo, agora2]
otu_list_full = pd.concat(otu_list)

otu_list_full.to_csv("05.01.20_all_otus_with_agora_or_mambo_models.csv")

#quick and dirty wins the race. 