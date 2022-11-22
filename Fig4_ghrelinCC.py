#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 16:40:17 2022

@author: albertsmith
"""

import pyDR
import numpy as np
import matplotlib.pyplot as plt
from misc_functions import load_helices
from pyDR.misc.Averaging import avgDataObjs

proj=pyDR.Project('Projects/ghrelin_iRED')

if not(len(proj['opt_fit'])):
    proj.detect.r_auto(7)
    proj.fit()['proc'].opt2dist(rhoz_cleanup=True)
    proj.save()


proj['opt_fit'].modes2bonds()
avgDataObjs(proj['iREDbond'])

data=proj['iREDbond']['.+AvOb'][0]


#%% Plot the results in chimeraX

proj['opt_fit']['.+rk1_.+ghrelin'].modes2bonds()
data=proj['iREDbond'][0]



command_line=proj.chimera.command_line

rho_index=5

hlx=[np.arange(41,71),np.arange(78,106),np.arange(114,146),np.arange(158,182),
     np.arange(209,238),np.arange(253,290),np.arange(296,327),np.arange(327,340)]


sel=data.select
#Fix representative selection
sel._mdmode=True
repr_sel=[sel.molsys.select_atoms('(resid {0} and name N CA HN) or (resid {1} and name C CA O)'.format(res,res-1)) for res in sel.sel1.resids[:-1]]
repr_sel.append(sel.molsys.select_atoms('resname SERO and name C3 C4 C5 C6 C7 C8'))
sel.repr_sel=repr_sel

commands=['turn x -90','turn y 190','turn z 5','view','set bgColor white',
     'lighting full','~show /A','show #10/A@N,C,CA','show #10/A:SERO&~@H*,O']

#%% Residues of Ghrelin
for k in range(13):
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
        sel.chimera(x=1.25*np.abs(data.CCnorm[rho_index][k][index]),index=index,color=clr)
        command_line(['~show #{0}'.format(m+1),
                      'show #{0}:'.format(m+1)+','.join([str(res) for res in resids])+'@N,C,CA',
                      'show #{0}:{1}@N'.format(m+1,resids[-1]+1),
                      'show #{0}:{1}@C'.format(m+1,resids[0]-1)])
    command_line(commands)
    command_line('color #10:{0}@N,CA|#10:{1}@C,CA black'.format(sel.sel1[k].resid,sel.sel1[k].resid-1))
    
    proj.chimera.savefig('rho{0}_resid{1}.png'.format(rho_index,sel.sel1[k].resid),options='transparentBackground True')

#%% Octanyl group from Ghrelin    
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
command_line('color #10:SERO@C3,C4,C5,C6,C7,C8 black'.format(sel.sel1[k].resid,sel.sel1[k].resid-1))

proj.chimera.savefig('rho{0}_octanyl.png'.format(rho_index,sel.sel1[k].resid),options='transparentBackground True')
    
    
