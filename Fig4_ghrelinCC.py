#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 16:40:17 2022

@author: albertsmith
"""

import pyDR
import numpy as np
import matplotlib.pyplot as plt
from pyDR.misc.Averaging import avgDataObjs

"If first time running pyDR with chimeraX, uncomment the lines below, and provide chimeraX path"
# from pyDR.chimeraX.chimeraX_funs import set_chimera_path
# set_chimera_path() #Put your own path to ChimeraX here!!
"e.g. set_chimera_path('/Applications/ChimeraX-1.2.5.app/Contents/MacOS/ChimeraX')"

proj=pyDR.Project('Projects/ghrelin_iRED')  #Load the project

if not(len(proj['opt_fit'])):    #Make sure the fully-processed data is in the project
    proj.detect.r_auto(7)
    proj.fit()['proc'].opt2dist(rhoz_cleanup=True)
    proj.save()


proj['opt_fit'].modes2bonds()   #Convert mode detector analysis to bond detector analysis
avgDataObjs(proj['iREDbond'])   #Average over the 3 trajectories



#%% Plot the results in chimeraX
data=proj['iREDbond'][-1]     #Data object from which we're plotting

command_line=proj.chimera.command_line   #Function to pass commands to ChimeraX

rho_index=5   #Which detector to show

hlx=[np.arange(41,71),np.arange(78,106),np.arange(114,146),np.arange(158,182),
     np.arange(209,238),np.arange(253,290),np.arange(296,327),np.arange(327,340)]   #Helix definitions used for color-coding

sel=data.select
#Fix representative selection
# The representative selection defines which atoms we encode the cross-correlation onto in ChimeraX
sel._mdmode=True
repr_sel=[sel.molsys.select_atoms('(resid {0} and name N CA HN) or (resid {1} and name C CA O)'.format(res,res-1)) for res in sel.sel1.resids[:-1]]
repr_sel.append(sel.molsys.select_atoms('resname SERO and name C3 C4 C5 C6 C7 C8'))
sel.repr_sel=repr_sel

commands=['turn x -90','turn y 190','turn z 5','view','set bgColor white',
     'lighting full','~show /A','show #10/A@N,C,CA','show #10/A:SERO&~@H*,O']

#%% Residues of Ghrelin
for k in range(13):      #Sweep over the 13 H–N containing residues in Ghrelin
    proj.chimera.close()
    for m in range(10):
        if m==9:   #if/elif/else state selects which atoms and which colors we show in each step
            resids=np.arange(2,16)
            clr=plt.get_cmap('Set1')(0)
        elif m==0:
            resids=np.concatenate(hlx)
            index=np.logical_not(np.isin(np.arange(33,340),resids))
            resids=np.arange(33,340)[index]
            clr=plt.get_cmap('tab10')(8)
        else:
            resids=hlx[m-1]
            clr=plt.get_cmap('tab10')(7 if m==8 else 7-m)
        index=np.argwhere(np.isin(sel.sel1.resids,resids)).squeeze()
        sel.chimera(x=1.25*np.abs(data.CCnorm[rho_index][k][index]),index=index,color=clr)   #Here we create the plot in ChimeraX
        command_line(['~show #{0}'.format(m+1),
                      'show #{0}:'.format(m+1)+','.join([str(res) for res in resids])+'@N,C,CA',
                      'show #{0}:{1}@N'.format(m+1,resids[-1]+1),
                      'show #{0}:{1}@C'.format(m+1,resids[0]-1)])
    command_line(commands)   #Change the chimeraX settings
    command_line('color #10:{0}@N,CA|#10:{1}@C,CA black'.format(sel.sel1[k].resid,sel.sel1[k].resid-1))  #Color the selected Ghrelin residue black
    
    proj.chimera.savefig(f'rho{rho_index}_resid{sel.sel1[k].resid}.png',options='transparentBackground True')  #Save the result



#%% Octanyl group from Ghrelin    
#Same as above, except for the octanyl group instead of the ghrelin backbone
proj.chimera.close()
for m in range(10):
    if m==9:
        resids=np.arange(2,16)
        clr=plt.get_cmap('Set1')(0)
    elif m==0:
        resids=np.concatenate(hlx)
        index=np.logical_not(np.isin(np.arange(33,340),resids))
        resids=np.arange(33,340)[index]
        clr=plt.get_cmap('tab10')(8)
    else:
        resids=hlx[m-1]
        clr=plt.get_cmap('tab10')(7 if m==8 else 7-m)
    index=np.argwhere(np.isin(sel.sel1.resids,resids)).squeeze()
    sel.chimera(x=1.25*np.abs(data.CCnorm[rho_index][-1][index]),index=index,color=clr)
    command_line(['~show #{0}'.format(m+1),
                  'show #{0}:'.format(m+1)+','.join([str(res) for res in resids])+'@N,C,CA',
                  'show #{0}:{1}@N'.format(m+1,resids[-1]+1),
                  'show #{0}:{1}@C'.format(m+1,resids[0]-1)])
command_line(commands)
command_line('color #10:SERO@C3,C4,C5,C6,C7,C8 black')

proj.chimera.savefig(f'rho{rho_index}_octanyl.png',options='transparentBackground True')
    

#%% Residues of Ghrelin: total correlation (not shown in paper)
k=4  #Just one residue as example. Could also be run in a for loop
#The results of this calculation show us why using the detector timescale filter is important for acquiring cross-correlation

proj.chimera.close()
for m in range(10):
    if m==9:
        resids=np.arange(2,16)
        clr=plt.get_cmap('Set1')(0)
    elif m==0:
        resids=np.concatenate(hlx)
        index=np.logical_not(np.isin(np.arange(33,340),resids))
        resids=np.arange(33,340)[index]
        clr=plt.get_cmap('tab10')(8)
    else:
        resids=hlx[m-1]
        clr=plt.get_cmap('tab10')(7 if m==8 else 7-m)
    index=np.argwhere(np.isin(sel.sel1.resids,resids)).squeeze()
    sel.chimera(x=1.25*np.abs(data.totalCCnorm[k][index]),index=index,color=clr)
    command_line(['~show #{0}'.format(m+1),
                  'show #{0}:'.format(m+1)+','.join([str(res) for res in resids])+'@N,C,CA',
                  'show #{0}:{1}@N'.format(m+1,resids[-1]+1),
                  'show #{0}:{1}@C'.format(m+1,resids[0]-1)])
command_line(commands)
command_line('color #10:{0}@N,CA|#10:{1}@C,CA black'.format(sel.sel1[k].resid,sel.sel1[k].resid-1))

proj.chimera.savefig(f'total_resid{sel.sel1[k].resid}.png',options='transparentBackground True')

plt.show()
    

#%% Save results to a text file (attached to paper, also in github in source_data folder)
from misc_functions import save_detectors

titles=['H–N correlation for bound GHSR']
data.label=np.array([lbl[:-5] for lbl in data.label],dtype=data.label.dtype)
data.label[-1]='octanyl'
save_detectors(fignum=4,data=[data],titles=titles,CCindex=[lbl for lbl in [*data.label[:13],data.label[-1]]])
