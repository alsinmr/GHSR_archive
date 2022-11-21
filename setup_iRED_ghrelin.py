#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 10:34:51 2022

@author: albertsmith
"""

import sys
sys.path.append('/work/home/alsi/python')
import pyDR
import os
import numpy as np
import matplotlib.pyplot as plt

def load_helices():
    name=list()
    resids=list()
    with open('/Users/albertsmith/Documents/Dynamics/GHSR/helix_assign.txt','r') as f:
        for line in f:
            name.append(line.strip().split(':')[0])
            ar=np.arange(\
                int(line.split('resid')[1].split('to')[0].strip()),
                int(line.split('to')[1].strip())+1)
            if len(resids):
                ar=ar[np.logical_not(np.isin(ar,np.concatenate(resids)))]
            resids.append(ar)
    return name,resids


#%% File locations
proj=pyDR.Project('Projects/ghrelin_iRED',create=True)


mddir='/work/public/ghrelin-receptor'
# mddir='/Volumes/My Book/GHSR'   

topo='WT-ghrelin_run1_just_protein.pdb'
trajs=[f'WT-ghrelin_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)]

#%% Loop over MD trajectories

for traj in trajs:
    molsys=pyDR.MolSys(topo=os.path.join(mddir,topo),
                       traj_files=os.path.join(mddir,traj),step=1,tf=355000,project=proj)
    
    sel=pyDR.MolSelect(molsys)
    sel.select_bond(Nuc='N')
    sel0=sel.molsys.select_atoms('resname SERO and name C3 C8')
    sel._mdmode=True
    sel.sel1+=sel0[0]
    sel.sel2+=sel0[1]
    
    frames=pyDR.Frames.FrameObj(sel)
    "Tensor frame"
    frames.tensor_frame(sel1=1,sel2=2)
    "Peptide plane frames"
    _,hlx=load_helices()
    "Overall"
    fr_sel=molsys.select_filter(resids=np.concatenate(hlx),filter_str='name CA')
    frames.new_frame(Type='superimpose',sel=fr_sel)
    
    
    data=frames.frames2iRED(rank=1)[1].iRED2data()
    data.detect.r_no_opt(15)
    data.fit(bounds=False)
    proj.save()
    

        
        


            
        






