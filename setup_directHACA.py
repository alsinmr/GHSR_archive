#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 10:24:26 2022

@author: albertsmith
"""

import sys
sys.path.append('/work/home/alsi/GitHub')
import pyDR
import os
import numpy as np
from pyDR.misc.Averaging import appendDataObjs
from misc_functions import load_helices



"""
When this is ported to the server, we need to change
1) The project location
2) The mddir
3) The time step
"""

#%% File locations
proj=pyDR.Project('Projects/backboneHACA',create=True)
mddir='/work/public/ghrelin-receptor'
# mddir='/Volumes/My Book/GHSR'

topos=['WT-apo_run1_just_protein.pdb','WT-ghrelin_run1_just_protein.pdb']
trajs=[[f'WT-apo_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)],
       [f'WT-ghrelin_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)]]


for topo,traj1 in zip(topos,trajs):
    for traj in traj1:
        sel=pyDR.MolSelect(topo=os.path.join(mddir,topo),
                           traj_files=os.path.join(mddir,traj),
                           step=100,project=proj,tf=355000)
        resids=sel.uni.residues.resids
        

        sel.select_bond(Nuc='CA',resids=[186,258,280])            #Select HA–CA bonds
        frames=pyDR.Frames.FrameObj(sel)    #Create frame object 
        frames.tensor_frame(sel1=1,sel2=2)  #Define tensor frame (HA–CA bond)
        _,hlx=load_helices()  #Use helices for alignment
        
        fr_sel=sel.molsys.select_filter(resids=np.concatenate(hlx),filter_str='name CA')
        frames.new_frame(Type='superimpose',sel=fr_sel)  #Remove overall motion
        
        frames.frames2data()
        
        proj.remove_data([-4,-3,-1],delete=True)

        proj[-1].detect.r_no_opt(15)
        proj[-1].fit()

        proj[-1].source.additional_info=None
        
        proj.update_info()
        proj.save()
            
