#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 10:26:48 2020

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

###Thanks to the excellent work of Building_microbiome_communities_03.10.20, we now have each OTU in Burns_2015 assigned a 
###model or model list. This code is to get a list of all the models I'll need to download to run these things. Then, I will 
#combine the models to make each OTU's 'composite model'

#take a look at species1/genus1 vs species2/genus 2: which one is more accurate/has best coverage?