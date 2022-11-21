#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 10:24:26 2022

@author: albertsmith
"""

import sys
sys.path.append('/work/home/alsi/python')
import pyDR
import os
import numpy as np
from pyDR.misc.Averaging import appendDataObjs


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

"""
When this is ported to the server, we need to change
1) The project location
2) The mddir
3) The time step
"""

#%% File locations
proj=pyDR.Project('Projects/backboneHN',create=True)
# proj=pyDR.Project('Projects/dummy',create=True)
mddir='/work/public/ghrelin-receptor'
# mddir='/Volumes/My Book/GHSR'

topos=['WT-apo_run1_just_protein.pdb','WT-ghrelin_run1_just_protein.pdb']
trajs=[[f'WT-apo_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)],
       [f'WT-ghrelin_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)]]


for topo,traj1 in zip(topos,trajs):
    for traj in traj1:
        sel=pyDR.MolSelect(topo=os.path.join(mddir,topo),
                           traj_files=os.path.join(mddir,traj),
                           step=1,project=proj,tf=355000)
        
        resids=sel.uni.residues.resids
        
        for k in range(5):  #Process in 5 chunks
            sel.select_bond(Nuc='N',resids=resids[62*k:62*(k+1)])            #Select H–N bonds
            frames=pyDR.Frames.FrameObj(sel)    #Create frame object 
            frames.tensor_frame(sel1=1,sel2=2)  #Define tensor frame (H–N bond)
            _,hlx=load_helices()  #Use helices for alignment
            
            fr_sel=sel.molsys.select_filter(resids=np.concatenate(hlx),filter_str='name CA')
            frames.new_frame(Type='superimpose',sel=fr_sel)  #Remove overall motion
            
            frames.frames2data()
            
            proj.remove_data([-4,-3,-1],delete=True)
            proj[-1].source.additional_info=f'chunk{k}'
            proj.update_info()
            proj[-1].detect.r_no_opt(15)
            proj[-1].fit()
            proj.clear_memory()
        
        appendDataObjs(proj[-5:])    
        proj.remove_data([-6,-5,-4,-3,-2],delete=True)
        proj[-1].source.additional_info=None
        
            
