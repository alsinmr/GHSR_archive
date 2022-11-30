#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 12:44:26 2022

@author: albertsmith
"""

import os
import matplotlib.pyplot as plt
import numpy as np

"""
These are results of chemical shift prediction and ring-current calculation from
SPARTA+. This results were ultimately not used in the paper, but we provide them
here for completeness

Reference:
    
SPARTA+: a modest improvement in empirical NMR chemical shift prediction by means of an artificial neural network
Yang Shen and Ad Bax, J. Biomol. NMR, 48, 13-22 (2010)
"""

#%% Make Ring current shift plot

rc_apo=list()
rc_ghrelin=list()
for file,rc in zip(['rc_apo','rc_ghrelin'],[rc_apo,rc_ghrelin]):
    with open(os.path.join('SPARTA',file+'.txt'),'r') as f:
        resids=np.array([int(res) for res in f.readline().strip().split()])
        for line in f:
            rc.append([float(x) for x in line.strip().split()])

rc_apo,rc_ghrelin=[np.array(rc) for rc in [rc_apo,rc_ghrelin]]

w=0.4
ax=plt.figure().add_subplot(111)
ax.bar(np.arange(len(resids))-w/2,rc_apo.mean(0),width=w,color='grey',edgecolor='black')
ax.bar(np.arange(len(resids))+w/2,rc_ghrelin.mean(0),width=w,color='darkgrey',edgecolor='black')
ax.set_xticks(range(len(resids)))
ax.set_xticklabels(resids)
ax.set_ylabel(r'Ring current $\delta(C\alpha)$')
ax.set_xlabel('Residue')

print('Resids (sorted from largest to smallest change in ring current effect)')
print(resids[np.argsort(np.abs(rc_apo.mean(0)-rc_ghrelin.mean(0)))][::-1])

#%% Make CS plot (not shown in text)

cs_apo=list()
cs_ghrelin=list()
for file,cs in zip(['shift_apo','shift_ghrelin'],[cs_apo,cs_ghrelin]):
    with open(os.path.join('SPARTA',file+'.txt'),'r') as f:
        resids=np.array([int(res) for res in f.readline().strip().split()])
        for line in f:
            cs.append([float(x) for x in line.strip().split()])

cs_apo,cs_ghrelin=[np.array(cs) for cs in [cs_apo,cs_ghrelin]]

w=0.4
ax=plt.figure().add_subplot(111)
ax.bar(np.arange(len(resids))-w/2,cs_apo.mean(0),width=w,color='grey',edgecolor='black')
ax.bar(np.arange(len(resids))+w/2,cs_ghrelin.mean(0),width=w,color='darkgrey',edgecolor='black')
ax.set_xticks(range(len(resids)))
ax.set_xticklabels(resids)
ax.set_ylabel(r'Pred. Chemical Shift $\delta(C\alpha)$')
ax.set_xlabel('Residue')
ax.set_ylim([58,63])
print('Resids (sorted from largest to smallest change in overall chemical shift')
print(resids[np.argsort(np.abs(cs_apo.mean(0)-cs_ghrelin.mean(0)))][::-1])