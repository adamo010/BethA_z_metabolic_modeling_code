#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 15:16:45 2020

@author: adamo010
"""

#goal: generate a txt file that contains a list of all pickle file names for MICOM communities from Burns 2015 data
import scipy
import numpy as np
import csv
import subprocess
import os
import pandas as pd
import fnmatch
import shutil
import glob
import re
import copy

os.chdir("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Pickle_files_for_comms")

filenames = glob.glob('*.pickle') #generates a list of file names
filenames

with open('pickle_file_list.txt', 'w') as f:
    for file in filenames:
        f.write("%s\n" % file)