#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 14:18:54 2023

@author: albertsmith
"""


import MDAnalysis as mda
import os
import numpy as np


#%% File Locations
mddir='/work/public/ghrelin-receptor'
mddir='/Volumes/My Book/GHSR/'

topos=['WT-apo_run1_just_protein.pdb','WT-ghrelin_run1_just_protein.pdb']
trajs=[[f'WT-apo_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)],
       [f'WT-ghrelin_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)]]

#%% Acquire reference distances (use the active, i.e. bound state as reference)

resids=[128,215,220,221,222,226,272,279,280,284,286,290,309,312]

r0=None
counter=0

topo=os.path.join(mddir,topos[1])
for traj in trajs[1][:1]:
    uni=mda.Universe(topo,os.path.join(mddir,traj))
    sel=uni.select_atoms('resid '+' '.join([str(res) for res in resids])+' and not name H* CA C N O')

    if r0 is None:
        r0=np.zeros(sel.positions.shape)

    for _ in uni.trajectory:
        counter+=1
        r0+=sel.positions
        
r0/=counter