#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 16:39:18 2022

@author: albertsmith
"""

from subprocess import run
import os
import MDAnalysis as mda


user='asmith'
host='calcium.nmrbox.org'

#%% Create pdbs
topos=['WT-apo_run1_just_protein.pdb','WT-ghrelin_run1_just_protein.pdb']
trajs=[[f'WT-apo_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)],
       [f'WT-ghrelin_run{k}_0.1ns_just_protein.xtc' for k in range(1,4)]]

mddir='/work/public/ghrelin-receptor'
working_dir='SHIFTX2'
if not(os.path.exists('SHIFTX2')):os.mkdir('SHIFTX2')

for topo,traj1 in zip(topos,trajs):
    for traj in traj1:
        print(traj)
        uni=mda.Universe(os.path.join(mddir,topo),os.path.join(mddir,traj))
        
        sel=uni.atoms if 'apo' in topo else uni.atoms.select_atoms('segid B')
        i=sel.residues.resnames=='HSD'
        sel.residues[i].resnames='HIS'  #Required for shiftX2
        
        title=traj.split('-')[1].rsplit('_',maxsplit=3)[0]+'_{}.pdb'
        
        for _ in uni.trajectory[:355000:500]: #Sample every 50 ns
            sel.write(os.path.join(working_dir,title.format(uni.trajectory.frame)))
        

#%% Run commands on NMRbox
cmd='mkdir /home/nmrbox/asmith/Desktop/shiftx_pdbs'
run(f'echo "{cmd}" | ssh {user}@{host}',shell=True)   #Create directory remotely

cmd=f"""scp {working_dir}/*.pdb {user}@{host}:/home/nmrbox/asmith/Desktop/shiftx_pdbs"""
run(cmd,shell=True)   #Copy the pdbs there

cmd=f"""echo "python2 /usr/software/SHIFTX2/shiftx2.py -b /home/nmrbox/asmith/Desktop/shiftx_pdbs/'*'.pdb" | ssh {user}@{host}"""
run(cmd,shell=True)  #Run shiftx on all pdbs

os.mkdir('cs')
cmd=f"""scp {user}@{host}:/home/nmrbox/asmith/Desktop/shiftx_pdbs/*.cs  {working_dir}/"""
run(cmd,shell=True)


#%% Collect the results
def collect(title='apo_run1'):
    shifts={'CA':{},'CB':{},'C':{}}

    for file in os.listdir(working_dir):
        if title in file:
            with open(os.path.join(working_dir,file),'r') as f:
                for line in f:
                    resid,resname,atomname,shift=line.strip.split(',')
                    if atomname in shifts:
                        pass