#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 16:18:57 2022

@author: albertsmith
"""

import sys
sys.path.append('/work/home/alsi/GitHub')
import pyDR
from pyDR.PCA import PCA
import numpy as np
from misc_functions import load_helices
import os
import matplotlib.py as plt


#%% Set up pca
mddir='/work/public/ghrelin-receptor'
# mddir='/Volumes/My Book/GHSR'

topos=['WT-apo_run1_just_protein.pdb','WT-ghrelin_run1_just_protein.pdb']
trajs=[[f'WT-apo_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)],
       [f'WT-ghrelin_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)]]

states=['apo','bound']

#%% Save part of PCA
""" 
We do the large processing (350000 time steps) on a server, but create the 
plots on a desktop, due to chimeraX requiring a window to run. Then, we rebuild 
the PCA objects based on saved files. 

Currently, PCA objects in pyDR cannot be saved.                                           
"""


if not(os.path.exists('PCA_results')):os.mkdir('PCA_results')

for topo,traj1,state in zip(topos,trajs,states):
    pca=PCA(pyDR.MolSelect(topo=os.path.join(mddir,topo),\
              traj_files=[os.path.join(mddir,traj) for traj in traj1[:1]],step=100)).\
              select_atoms('name N C CA')
    with open(os.path.join('PCA_results',f'{state}_covar.dat'),'wb') as f:
        np.save(f,pca.CoVar,allow_pickle=False)
    with open(os.path.join('PCA_results',f'{state}_pcamp.data'),'wb') as f:
        np.save(f,pca.PCamp,allow_pickle=False)
    with open(os.path.join('PCA_results',f'{state}_mean.data'),'wb') as f:
        np.save(f,pca.mean,allow_pickle=False)
