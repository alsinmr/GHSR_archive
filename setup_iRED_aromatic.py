#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 17:53:20 2022

@author: albertsmith
"""

import sys
sys.path.append('/work/home/alsi/GitHub')
import pyDR
import os
import numpy as np
from misc_functions import helix_only


mddir='/work/public/ghrelin-receptor'
# mddir='/Volumes/My Book/GHSR'

topos=['WT-apo_run1_just_protein.pdb','WT-ghrelin_run1_just_protein.pdb']
trajs=[[f'WT-apo_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)],
       [f'WT-ghrelin_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)]]

proj=pyDR.Project('Projects/aromatic_iRED',create=True)

resids=helix_only()[::3]

for topo,traj1 in zip(topos,trajs):
    for traj in traj1:
        select=pyDR.MolSelect(topo=os.path.join(mddir,topo),
        traj_files=os.path.join(mddir,traj),step=1,project=proj,tf=355000)
        
        segid='A B' if 'ghrelin' in topo else 'A'
        select_atoms=select.molsys.select_atoms
        sel1=select_atoms('resname PHE TYR and name CZ and segid {0}'.format(segid))+\
            select_atoms('resname HSD and name NE2 and segid {0}'.format(segid))+\
                select_atoms('resname TRP and name CZ3 and segid {0}'.format(segid))
        sel2=select_atoms('resname PHE TYR HSD and name CG and segid {0}'.format(segid))+\
            select_atoms('resname TRP and name CD1 and segid {0}'.format(segid))
        sel3=select_atoms('resname PHE TYR and name CD1 and segid {0}'.format(segid))+\
            select_atoms('resname HSD and name ND1 and segid {0}'.format(segid))+\
                select_atoms('resname TRP and name NE1 and segid {0}'.format(segid))
                
        i=sel1.resids.argsort()
        
        sel1,sel2,sel3=sel1[i],sel2[i],sel3[i]
        
        
        
        sel1,sel2,sel3=[sel+select_atoms('name '+name+' and resid '+ ' '.join([str(r) for r in resids])) \
                                         for sel,name in zip([sel1,sel2,sel3],['N','CA','C'])]
        
                
        sel0=select_atoms('name CA and resid '+' '.join([str(r) for r in resids]))
                
        repr_sel=[s.residue.atoms for s in sel1]
        
        select.sel1,select.sel2,select.repr_sel=[sel1,sel2,repr_sel]
        select.label=np.array([lbl.split('_')[0]+(' BB' if 'N_CA' in lbl else '') for lbl in select.label],dtype=str)
        
        frames=pyDR.Frames.FrameObj(select)
        frames.tensor_frame(Type='bond',sel1=sel1,sel2=sel2,sel3=sel3)
        
        ired=frames.md2iRED(rank=1)
        ired.iRED2data()
        
        proj[-1].detect.r_no_opt(15)
        proj[-1].fit()
        proj.save()
