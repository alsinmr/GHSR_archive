#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 15:09:47 2023

@author: albertsmith
"""

import pyDR
from pyDR.PCA import PCA
import numpy as np
import os
import matplotlib.pyplot as plt

"If first time running pyDR with chimeraX, uncomment the lines below, and provide chimeraX path"
# from pyDR.chimeraX.chimeraX_funs import set_chimera_path
# set_chimera_path() #Put your own path to ChimeraX here!!
"e.g. set_chimera_path('/Applications/ChimeraX-1.2.5.app/Contents/MacOS/ChimeraX')"

#%% File locations

topos=['WT-apo_run1_just_protein.pdb','WT-ghrelin_run1_just_protein.pdb']
states=['apo','bound']

#%% Plot path for apo and bound form
fig,ax=plt.subplots(2,3)
cmap=plt.get_cmap('jet')
lengths=[[372685,364134,328181],[364893,357953,342154]]
step=100

for topo,state,ax0,len0 in zip(topos,states,ax,lengths):
    #Re-load the PCA data
    pca=PCA(pyDR.MolSelect(topo=topo)).select_atoms('name N C CA')
    with open(os.path.join('PCA_results',f'{state}_covar.data'),'rb') as f:
        pca._covar=np.load(f,allow_pickle=False)
    with open(os.path.join('PCA_results',f'{state}_pcamp.data'),'rb') as f:
        pca._pcamp=np.load(f,allow_pickle=False)
    with open(os.path.join('PCA_results',f'{state}_mean.data'),'rb') as f:
        pca._mean=np.load(f,allow_pickle=False)    
    with open(os.path.join('PCA_results',f'{state}_Lambda.data'),'rb') as f:
        pca._lambda=np.load(f,allow_pickle=False) 
    with open(os.path.join('PCA_results',f'{state}_PC.data'),'rb') as f:
        pca._PC=np.load(f,allow_pickle=False)
        
        
    for k,a in enumerate(ax0):
        i=range(sum(len0[:k]),sum(len0[:k+1]),step)
        a.scatter(pca.PCamp[0,i],pca.PCamp[1,i],3,c=cmap(np.linspace(0,1,len(i))))
        a.set_xlabel('PC 0')
        a.set_ylabel('PC 1')
        a.set_xlim([-150,150])
        a.set_ylim([-150,150])
        a.text(-100,100,f'{state}\nrun {k+1}',verticalalignment='top')
fig.set_size_inches([10.5,7])
fig.tight_layout()        