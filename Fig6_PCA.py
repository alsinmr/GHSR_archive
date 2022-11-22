#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 15:00:31 2022

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


server=True    #Use this as a flag to decide if we're on the server or not

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

if server:   #Delete or set to True if you actually want to run these lines
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



#%% Re-create the PCA object from file

if not(server):
    #String for selections in chimeraX
    hlx_sel='setattr :'+','.join(str(int(i)) for i in \
            np.concatenate([res if 'in' in name or 'ex' in name or name=='8' else [] for name,res in zip(*load_helices())]))+\
            ' res ss_type 1'
            
    cmx_commands=['lighting soft','graphics silhouettes true',hlx_sel,'ribbon']
    
    for topo,state in zip(topos,states):
        pca=PCA(pyDR.MolSelect(topo=os.path.join(mddir,topo))).select_atoms('name N C CA')
        with open(os.path.join('PCA_results',f'{state}_covar.dat'),'rb') as f:
            pca._covar=np.load(f,allow_pickle=False)
        with open(os.path.join('PCA_results',f'{state}_pcamp.data'),'rb') as f:
            pca._pcamp=np.load(f,allow_pickle=False)
        with open(os.path.join('PCA_results',f'{state}_mean.data'),'rb') as f:
            pca._mean=np.load(f,allow_pickle=False)
            
        "Plot the 1st and 2nd components of the PCA"
        fig=plt.figure()
        ax=fig.add_subplot(111)
        hdl=pca.plot(n0=0,n1=1,maxbin=120,cmap='nipy_spectral',nbins=100,ax=ax)
        ax.set_xlabel('PC 1')
        ax.set_ylabel('PC 2')
        hdl1=plt.colorbar(hdl[-1],ax=ax)
        hdl1.set_label('count')
        fig.set_size_inches([5.8,4.4])
        fig.tight_layout()
        fig.savefig(f'PCA_results/PC1_2_{state}.pdf')       
        
        "Plot in chimera"
        for k in range(2):
            pca.project.chimera.current=k
            pca.chimera(n=k,std=1)
            pca.project.chimera.command_line(cmx_commands)

        "Make free energy plots"
        fig=plt.figure()
        ax=fig.add_subplot(111)
        
        fig.set_size_inches([5.8,4.4])
        fig.tight_layout()
        fig.savefig(f'PCA_results/free_energy_{state}.pdf')
        
        mb=120
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
        ax.set_ylim([0,25])
        
        

