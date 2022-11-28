#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 16:27:42 2022

@author: albertsmith
"""

import pyDR
from misc_functions import load_helices
from pyDR.misc.Averaging import avgDataObjs


#%% Finish processing
proj=pyDR.Project('Projects/CACB')

#We'll check if the operations have already been run/saved before executing them

#Average the 3x2 trajectories
if not(len(proj['.+AvOb'])):
    avgDataObjs(proj['.+apo'])
    avgDataObjs(proj['.+ghrelin'])
    
    
if not(len(proj['opt_fit'])):  #Convert no-opt detectors to 7 optimized detectors
    proj['no_opt']['.+AvOb'].detect.r_auto(7)
    proj['no_opt']['.+AvOb'].fit()
    proj['proc'].opt2dist(rhoz_cleanup=True)
    proj.save()

#%% Part A: Zoom in on Loop 

files=['WT-apo_run1_0.1ns_just_protein','WT-ghrelin_run1_0.1ns_just_protein']

sub=proj['opt_fit']

_,hlx=load_helices()
ranges=[[180,200],[247,328]]
scaling=[2,4]
# index=[[np.argwhere(np.logical_and(s.label>=r[0],s.label<=r[1]))[:,0] for s in sub] for r in ranges]

proj.chimera.saved_commands=['set bgColor white','~show','show @N,C,CA,CB',
                             'lighting full','graphics silhouettes false','clip off',
                             'view']

r=ranges[0]
sc=scaling[0]
for k in range(7):
    proj.chimera.close()
    proj['opt_fit'].chimera(rho_index=k,scaling=sc)
    
    proj.chimera.command_line(['align #2/B@CA toAtoms #1@CA','align #4/B@CA toAtoms #3@CA',
                               'view #1:180-200|#1:288-297',
                               'zoom .9','move y -5','~show #2|#4','clip off','clip far 10'])
    
    proj.chimera.savefig('HIS186_apo_rho{0}.png'.format(k),options='transparentBackground True')
    
    proj.chimera.command_line(['~show #1','show #2@N,C,CA,CB','clip off','clip far 10','color #2/A slate grey'])
    proj.chimera.savefig('HIS186_bound_rho{0}.png'.format(k),options='transparentBackground True')

    
#%% Part B: Zoom in on helices 4-7

r=ranges[1]
sc=2
r0=[207,256]
rho_index=5

# ranges=[[157,184],[210,239],[253,327]]
ranges=[[201,239],[253,327]]

sel_str=':'+','.join(['{}-{}'.format(*r) for r in ranges])

proj.chimera.close()
proj['opt_fit'].chimera(rho_index=rho_index,scaling=sc)

proj.chimera.command_line(['~show','~ribbon','sel #1'+sel_str+'@N,C,CA,CB',
                           'ribbon #1','~ribbon sel','show sel',
                           'color #1&~sel light steel blue','~sel',
                           'turn x -90','turn y -75','turn z 15','view']) 
proj.chimera.savefig(f'Hlx5to7_apo_rho{rho_index}.png',options='transparentBackground True')

proj.chimera.command_line(['~show','~ribbon','sel #2'+sel_str+'@N,C,CA,CB',
                           'ribbon #2','~ribbon sel','show sel',
                           'color #2&~sel light steel blue','~sel','view #2','clip off',
                           'zoom 1.3','turn y -10','color #2/A slate grey'])
proj.chimera.savefig(f'Hlx5to7_bound_rho{rho_index}.png',options='transparentBackground True')
    