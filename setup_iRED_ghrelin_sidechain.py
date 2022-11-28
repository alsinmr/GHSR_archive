#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 10:34:51 2022

@author: albertsmith
"""

import sys
sys.path.append('/work/home/alsi/GitHub')
import pyDR
import os
import numpy as np
from misc_functions import load_helices


#%% File locations
proj=pyDR.Project('Projects/ghrelin_iRED_sidechain',create=True)


mddir='/work/public/ghrelin-receptor'
# mddir='/Volumes/My Book/GHSR'   

topo='WT-ghrelin_run1_just_protein.pdb'
trajs=[f'WT-ghrelin_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)]

#%% Loop over MD trajectories

for traj in trajs:
    molsys=pyDR.MolSys(topo=os.path.join(mddir,topo),
                       traj_files=os.path.join(mddir,traj),step=1,tf=355000,project=proj)
    
    sel=pyDR.MolSelect(molsys)
    sel.select_bond(Nuc='sidechain')
    sel0=sel.molsys.select_atoms('resname SERO and name C3 C8')
    sel._mdmode=True
    sel.sel1+=sel0[0]
    sel.sel2+=sel0[1]
    sel.repr_sel=[s.residue.atoms for s in sel.sel1]
    sel.repr_sel[-1]=sel.repr_sel[-1][18:]
    sel.repr_sel[2]=sel.repr_sel[2][:18]
    
    
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
    

        
        


            
        






