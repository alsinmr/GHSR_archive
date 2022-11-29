#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 15:00:31 2022

@author: albertsmith
"""



import pyDR
from pyDR.PCA import PCA
import numpy as np
from misc_functions import load_helices
import os
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh


#%% File locations
# mddir='/work/public/ghrelin-receptor'
mddir='/Volumes/My Book/GHSR'

topos=['WT-apo_run1_just_protein.pdb','WT-ghrelin_run1_just_protein.pdb']
states=['apo','bound']

#%% Re-create the PCA object from file

#String for selections in chimeraX
hlx_sel='setattr :'+','.join(str(int(i)) for i in \
        np.concatenate([res if 'in' in name or 'ex' in name or name=='8' else [] for name,res in zip(*load_helices())]))+\
        ' res ss_type 1'
        
cmx_commands=['lighting soft','graphics silhouettes true',hlx_sel,'ribbon']

for topo,state in zip(topos,states):
    """
    Note that the pyDR PCA module is not fully set up for export, so this is
    a little hack-y. We have saved parts of the analysis on a server and
    reload them into the hidden variables in pca to restore functionality.
    """
    pca=PCA(pyDR.MolSelect(topo=os.path.join(mddir,topo))).select_atoms('name N C CA')
    with open(os.path.join('PCA_results',f'{state}_covar.dat'),'rb') as f:
        pca._covar=np.load(f,allow_pickle=False)
    with open(os.path.join('PCA_results',f'{state}_pcamp.data'),'rb') as f:
        pca._pcamp=np.load(f,allow_pickle=False)
    with open(os.path.join('PCA_results',f'{state}_mean.data'),'rb') as f:
        pca._mean=np.load(f,allow_pickle=False)    
    w,v=eigsh(pca.CoVar,k=10,which='LM')    
    i=np.argsort(w)[::-1]
    pca._lambda,pca._PC=w[i],v[:,i]

    "Plot the 1st and 2nd components of the PCA"
    fig=plt.figure()
    ax=fig.add_subplot(111)
    hdl=pca.plot(n0=0,n1=1,maxbin=120,cmap='nipy_spectral',nbins=100,ax=ax)
    ax.set_xlabel('PC 0')
    ax.set_ylabel('PC 1')
    hdl1=plt.colorbar(hdl[-1],ax=ax)
    hdl1.set_label('count')
    fig.set_size_inches([5.8,4.4])
    fig.tight_layout()
    fig.savefig(f'PCA_results/PC1_2_{state}.pdf')       
    
    
    pca.project.chimera.saved_commands=cmx_commands
    pca.hist2struct(nmax=3,maxbin=120,cmap='nipy_spectral',nbins=100)
    
    # "Draw in chimera"
    # chimera=pca.project.chimera
    # for k in range(2):
    #     chimera.current=k
    #     pca.chimera(n=k,std=1)
    #     chimera.command_line(cmx_commands)
    #     if state=='bound':
    #         chimera.command_line(['color #1/A blue','color #2/A yellow'])
            
    "Make 3D plots (not shown in paper)"
    fig=plt.figure()
    ax=[fig.add_subplot(1,3,k+1,projection='3d') for k in range(3)]
    step=pca.PCamp.shape[1]//3
    cmap=plt.get_cmap('jet')
    lengths=[372685,364133,362956] if state=='apo' else [364892,357952,367384]
    for k,a in enumerate(ax):
        xyz=pca.PCamp[:3,sum(lengths[:k]):sum(lengths[:k+1]):100]
        c=cmap(np.linspace(0,1,xyz.shape[1]))
        a.scatter3D(*xyz,s=5,c=c)
    for a in ax:
        a.set_xlim([-120,120]),a.set_ylim([-120,120]),a.set_zlim([-70,70])
        a.set_xlabel('PC 0'),a.set_ylabel('PC 1'),a.set_zlabel('PC 2')
    fig.set_size_inches([16.6,4.9])
    fig.tight_layout()

    "Make free energy plots (not shown in paper– not really accurate)"
    fig=plt.figure()
    ax=fig.add_subplot(111)
    
    mb=140
    bins=np.linspace(-mb,mb,200)
    x=(bins[1:]+bins[:-1])/2
    RT=8.314*298/1000
    
    A=np.histogram(pca.PCamp,bins=bins)[0].astype(float)
    A/=A.sum()
    DelG=-np.log(A/A.max())*RT
    ax.plot(x,DelG)
    ax.set_xlabel('PCA 1')
    ax.set_ylabel(r'$\Delta$G / (kJ/mol)')
    ax.set_xlim([-mb,mb])
    ax.set_ylim([0,35])
    
    fig.set_size_inches([5.8,4.4])
    fig.tight_layout()
    fig.savefig(f'PCA_results/free_energy_{state}.pdf')
    
        

