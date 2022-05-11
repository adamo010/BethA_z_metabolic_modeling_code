#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 11:11:15 2020

@author: adamo010
"""
#today I am going to work on building metabolic models from scratch using CarveMe.
#step 0: pip install carveme into this environment. 
#also need to install cplex; the default is (hopefully, according to update release notes) gurobi, but carveme won't
#import without cplex
#step 0.1: pip install cplex
#oh boy, have to get an academic license and everything for this sucker. Great. Stay tuned. 
#cplex is now installed (minus the API, which I can't figure out), as is CarveMe

import cplex
import carveme
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

carve (Victivallis_vadensis_DSM_14823_protein.faa)
#oh. What if this goes in the command line and not in spyder?

