#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 12:49:04 2023

@author: albertsmith
"""

import os
import numpy as np
import matplotlib.pyplot as plt






#%% Calculate Q (276 compared to aromatics)
with open(os.path.join('contact_score','R.data'),'rb') as f:
    R=np.load(f,allow_pickle=False)
with open(os.path.join('contact_score','r0.data'),'rb') as f:
    r0=np.load(f,allow_pickle=False)

t=np.arange(R.shape[0]//3)*10

beta=5

fig,ax=plt.subplots(3,2)
for Lambda,ax0 in zip([1.8,1],ax.T):
    
    Q=(1/(1+np.exp(beta*(R-Lambda*r0)))).mean(1)
    
    for k,a in enumerate(ax0):
        a.plot(t/1e3,Q[len(Q)//3*k:len(Q)//3*(k+1)])
        if a.get_subplotspec().is_last_row():
            a.set_xlabel(r't / $\mu$s')
        else:
            a.set_xticklabels('')
        if a.get_subplotspec().is_first_row():
            a.set_title(r'$\lambda = $'+f'{Lambda}')
        if a.get_subplotspec().is_first_col():
            a.set_ylabel('Q')
        a.text(30,a.get_ylim()[0]+np.diff(a.get_ylim())*.05,f'run {k+1}')
        
        
#%% Calculate Q (aromatics vs. aromatics)
with open(os.path.join('contact_score','R_all.data'),'rb') as f:
    R=np.load(f,allow_pickle=False)
with open(os.path.join('contact_score','r0_all.data'),'rb') as f:
    r0=np.load(f,allow_pickle=False)


t=np.arange(R.shape[0]//3)*10

beta=5

fig,ax=plt.subplots(3,2)
for Lambda,ax0 in zip([1.8,1],ax.T):
    
    Q=(1/(1+np.exp(beta*(R-Lambda*r0)))).mean(1)
    
    for k,a in enumerate(ax0):
        a.plot(t/1e3,Q[len(Q)//3*k:len(Q)//3*(k+1)])
        if a.get_subplotspec().is_last_row():
            a.set_xlabel(r't / $\mu$s')
        else:
            a.set_xticklabels('')
        if a.get_subplotspec().is_first_row():
            a.set_title(r'$\lambda = $'+f'{Lambda}')
        if a.get_subplotspec().is_first_col():
            a.set_ylabel('Q')
        a.text(30,a.get_ylim()[0]+np.diff(a.get_ylim())*.05,f'run {k+1}')