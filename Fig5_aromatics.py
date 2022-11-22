#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 18:58:35 2022

@author: albertsmith
"""

import pyDR
import numpy as np
import matplotlib.pyplot as plt
from pyDR.misc.Averaging import avgDataObjs
from copy import copy

proj=pyDR.Project('Projects/aromatic_iRED')

#%% FInish processing
if not(len(proj['opt_fit'])):
    proj.detect.r_auto(7)
    proj.fit()['proc'].opt2dist(rhoz_cleanup=True)
    proj.save()
    
proj['opt_fit'].modes2bonds()

avgDataObjs(proj['iREDbond']['.+apo'])
avgDataObjs(proj['iREDbond']['.+ghrelin'])

sub=proj['iREDbond']['rk1']['.+AvOb']
sub=proj['iREDbond']['rk1']


resids=[128,215,220,221,222,226,272,279,280,284,286,290,309,312]
sel=copy(sub[1].select)
sel.select_bond(Nuc='15N',resids=resids)
sel.repr_sel=[s.residues.atoms for s in sel.sel1]

index0=[np.argwhere(s.label=='276')[0,0] for s in sub]
index=[np.isin(s.label,[str(r) for r in resids]) for s in sub]  #These are the residues we'll correlate with

#%% Plot set up
nd=sub[0].ne-1
w=0.45
cmap=plt.get_cmap('tab10')

#%% Make bar plot with just 276
ax=plt.figure().add_subplot(111)
hatch=['','///']
for m,(i0,s,h) in enumerate(zip(index0,sub,hatch)):
    for k,R in enumerate(s.R[i0]):
        ax.bar(k-w/2+m*w,R,width=w,color=cmap(k),hatch=h)