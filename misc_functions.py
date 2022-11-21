#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:56:23 2022

@author: albertsmith
"""

import numpy as np


#%% Defines where helices are in GHSR
def load_helices():
    name=list()
    resids=list()
    with open('helix_assign.txt','r') as f:
        for line in f:
            name.append(line.strip().split(':')[0])
            ar=np.arange(\
                int(line.split('resid')[1].split('to')[0].strip()),
                int(line.split('to')[1].strip())+1)
            if len(resids):
                ar=ar[np.logical_not(np.isin(ar,np.concatenate(resids)))]
            resids.append(ar)
    return name,resids
