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


def helix_only(select=None):
    names,hlx=load_helices()
    i=np.arange(32,340)
    i0=np.concatenate(np.array(hlx)[np.array(['L' in name or 'NT' in name for name in names],dtype=bool)])
    i0=np.logical_not(np.isin(i,i0))
    i=i[i0]
    if select is None:return i
    return np.array([s.resids[0] in i for s in select.sel1],dtype=bool)