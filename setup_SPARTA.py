#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 09:03:42 2022

@author: albertsmith
"""

import MDAnalysis as mda
import os
import subprocess


#%% Load trajectories and save truncated trajectory
mddir='/work/public/ghrelin-receptor'
#mddir='/Volumes/My Book/GHSR'

topos=['WT-apo_run1_just_protein.pdb','WT-ghrelin_run1_just_protein.pdb']
trajs=[[f'WT-apo_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)],
       [f'WT-ghrelin_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)]]

if not(os.path.exists('SPARTA')):os.mkdir('SPARTA')

for topo,traj1 in zip(topos,trajs):
    save_name=os.path.join('SPARTA',topo.split('_')[0]+'_4sparta.xtc')
    
    if not(os.path.exists(save_name)):
        uni=mda.Universe(os.path.join(mddir,topo),*[os.path.join(mddir,traj) for traj in traj1[:1]])
        uni.residues[uni.residues.resnames=='HSD'].resnames='HIS'
        
        step=len(uni.trajectory)//2500
        
        with mda.Writer(save_name,uni.atoms.n_atoms,step=step) as W:
            for ts in uni.trajectory[::step]:
                W.write(uni.atoms)
                
                
                
                
#%% Now run SPARTA for the truncated trajectories (~2500 frames each)
resids=[str(res) for res in [128, 215, 220, 221, 222, 226, 272, 276, 279, 280, 284, 286, 290, 309, 312]]

for topo in topos:
    save_name=os.path.join('SPARTA',topo.split('_')[0]+'_4sparta.xtc')
    
    uni=mda.Universe(os.path.join(mddir,topo),save_name)
    
    name=topo.split('_')[0].split('-')[1]
    file_rc=os.path.join('SPARTA',f'rc_{name}.txt')
    file_shift=os.path.join('SPARTA',f'shift_{name}.txt')
    
    with open(file_rc,'w') as frc,open(file_shift,'w') as fs:
        frc.write('\t'.join([str(res) for res in resids])+'\n')
        fs.write('\t'.join([str(res) for res in resids])+'\n')
        
        for k,_ in enumerate(uni.trajectory):
            if k%10==0:print(k)
            sel=uni.select_atoms('protein')
            sel.residues[sel.residues.resnames=='HSD'].resnames='HIS'
            sel.write(os.path.join('SPARTA','structure.pdb'))
            subprocess.call(['/bin/csh','-c','sparta+ -in SPARTA/structure.pdb'],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            with open('pred.tab','r') as f:
                for line in f:
                    if len(line.strip().split())==9:
                        resid,_,atom,_,shift,_,ring,_,_=line.strip().split()
                        if resid in resids and atom=='CA':
                            fs.write(shift+'\t')
                            frc.write(ring+'\t')
            fs.write('\n')
            frc.write('\n')
os.remove('pred.tab')
os.remove('struct.tab')      
