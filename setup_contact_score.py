#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 10:28:49 2023

@author: albertsmith
"""

import numpy as np
import MDAnalysis as mda
import os
from copy import copy

mddir='/work/public/ghrelin-receptor'

topos=['WT-apo_run1_just_protein.pdb','WT-ghrelin_run1_just_protein.pdb']

trajs=[[f'WT-apo_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)],
       [f'WT-ghrelin_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)]]


# mddir='/Volumes/My Book/GHSR'
# topos=['WT-apo_run1_just_protein.pdb','WT-ghrelin_run1_just_protein.pdb'][::-1]
# trajs=[['WT-apo_run1_0.1ns_just_protein.xtc'],['WT-ghrelin_run1_0.1ns_just_protein.xtc']][::-1]


uni=mda.Universe(os.path.join(mddir,topos[1]),[os.path.join(mddir,traj0) for traj0 in trajs[1]])


#%% Prepare selections
sel1=uni.select_atoms('resid 276 and not name H*')


resids=[128,215,220,221,222,226,272,279,280,284,286,290,309,312]
sel2=[uni.select_atoms(f'resid {res} and not name H*') for res in resids]
index=[]
for k,sel20 in enumerate(sel2):
    start=index[-1][-1]+1 if len(index) else 0
    index.append(np.arange(start,start+len(sel20)))

sel2=sum(sel2)
sel1=copy(sel2)

l1=len(sel1)
sel1=sum([sel1 for _ in range(len(sel2))])
sel2=sum([sum([sel20 for _ in range(l1)]) for sel20 in sel2])

#%% Sweep over bound (active) trajectory
step=100
r=np.zeros([len(uni.trajectory)//step+1,len(sel1)])


for k,_ in enumerate(uni.trajectory[::step]):
    if k%1000==0:print(f'{k} out of {len(uni.trajectory[::step])}')
    r[k]=np.sqrt(((sel1.positions-sel2.positions)**2).sum(1))
        
#%% Determine mean distances, index of contacts
r0=r.mean(0)
i=np.logical_and(np.logical_and(r0<4.5,np.abs(sel2.resids-sel1.resids)>=3),r0>.1)

#%% Prepare selections (bound)


tf=355000
step=100
R=[]

for traj0 in trajs[1]:
    uni=mda.Universe(os.path.join(mddir,topos[1]),os.path.join(mddir,traj0))
    sel1=uni.select_atoms('resid 276 and not name H*')
    
    
    resids=[128,215,220,221,222,226,272,279,280,284,286,290,309,312]
    sel2=[uni.select_atoms(f'resid {res} and not name H*') for res in resids]
    index=[]
    for k,sel20 in enumerate(sel2):
        start=index[-1][-1]+1 if len(index) else 0
        index.append(np.arange(start,start+len(sel20)))
    
    sel2=sum(sel2)
    sel1=copy(sel2)
    
    l1=len(sel1)
    sel1=sum([sel1 for _ in range(len(sel2))])
    sel2=sum([sum([sel20 for _ in range(l1)]) for sel20 in sel2])
    
    sel1=sel1[i]
    sel2=sel2[i]
    
    #%% Sweep over the trajectory
    R.append(np.zeros([len(uni.trajectory[:tf:step]),len(sel1)]))
    for k,_ in enumerate(uni.trajectory[:tf:step]):
        if k%1000==0:print(f'{k} out of {len(uni.trajectory[:tf:step])}')
        R[-1][k]=np.sqrt(((sel1.positions-sel2.positions)**2).sum(1))
R=np.concatenate(R,axis=0)

#%% Save results
r0=r0[i]
if not(os.path.exists('contact_score')):os.mkdir('contact_score')
with open(os.path.join('contact_score','R.data'),'wb') as f:
    np.save(f,R,allow_pickle=False)
with open(os.path.join('contact_score','r0.data'),'wb') as f:
    np.save(f,r0,allow_pickle=False)



