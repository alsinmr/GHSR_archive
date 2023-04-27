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
import os


#%% Set up pca
mddir='/work/public/ghrelin-receptor'

topos=['WT-apo_run1_just_protein.pdb','WT-apo_run1_just_protein.pdb','WT-apo_run1_just_protein.pdb',
       'WT-ghrelin_run1_just_protein.pdb','WT-ghrelin_run1_just_protein.pdb','WT-ghrelin_run1_just_protein.pdb']

trajs=[*[f'WT-apo_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)],
       *[f'WT-ghrelin_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)]]


# mddir='/Volumes/My Book/GHSR'
# topos=['WT-apo_run1_just_protein.pdb','WT-ghrelin_run1_just_protein.pdb']
# trajs=['WT-apo_run1_0.1ns_just_protein.xtc','WT-ghrelin_run1_0.1ns_just_protein.xtc']

#%% Save part of PCA
""" 
We do the large processing (350000 time steps) on a server, but create the 
plots on a desktop, due to chimeraX requiring a window to run. Then, we rebuild 
the PCA objects based on saved files. 

Currently, PCA objects in pyDR cannot be saved.                                           
"""


if not(os.path.exists('PCA_avb_results')):os.mkdir('PCA_avb_results')

tf=355000
step=1
pos=None


for topo,traj in zip(topos,trajs):
    pca=PCA(pyDR.MolSelect(topo=os.path.join(mddir,topo),\
              traj_files=os.path.join(mddir,traj),step=step,tf=tf)).\
              select_atoms('name N C CA and resid 30-400')
              
    if pos is None:
        pos=pca.pos
    else:
        pos=np.concatenate((pos,pca.pos),axis=0)
pca._pos=pos
pca.align()

pca.runPCA(n=10)
with open(os.path.join('PCA_avb_results','covar.data'),'wb') as f:
    np.save(f,pca.CoVar,allow_pickle=False)
with open(os.path.join('PCA_avb_results','Lambda.data'),'wb') as f:
    np.save(f,pca.Lambda,allow_pickle=False)
with open(os.path.join('PCA_avb_results','PC.data'),'wb') as f:
    np.save(f,pca.PC,allow_pickle=False)
with open(os.path.join('PCA_avb_results','pcamp.data'),'wb') as f:
    np.save(f,pca.PCamp,allow_pickle=False)
with open(os.path.join('PCA_avb_results','mean.data'),'wb') as f:
    np.save(f,pca.mean,allow_pickle=False)
