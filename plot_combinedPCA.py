#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 13:53:43 2023

@author: albertsmith
"""

import pyDR
from pyDR.PCA import PCA
import numpy as np
import os

topo='WT-apo_run1_just_protein.pdb'

pca=PCA(pyDR.MolSelect(topo=topo)).select_atoms('name N C CA')

with open(os.path.join('PCA_avb_results','covar.data'),'rb') as f:
    pca._covar=np.load(f,allow_pickle=False)
with open(os.path.join('PCA_avb_results','pcamp.data'),'rb') as f:
    pca._pcamp=np.load(f,allow_pickle=False)
with open(os.path.join('PCA_avb_results','mean.data'),'rb') as f:
    pca._mean=np.load(f,allow_pickle=False)
with open(os.path.join('PCA_avb_results','Lambda.data'),'rb') as f:
    pca._lambda=np.load(f,allow_pickle=False) 
with open(os.path.join('PCA_avb_results','PC.data'),'rb') as f:
    pca._PC=np.load(f,allow_pickle=False)
pca._t=np.linspace(0,6*355000*.1,6*355000)

fig=pca.hist2struct(nmax=3,maxbin=150,cmap='nipy_spectral',nbins=100,ref_struct=False)[0].figure