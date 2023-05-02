#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 12:49:04 2023

@author: albertsmith
"""

import os
import numpy as np
import matplotlib.pyplot as plt


with open(os.path.join('contact_score','R.data'),'rb') as f:
    R=np.load(f,allow_pickle=False)
with open(os.path.join('contact_score','r0.data'),'rb') as f:
    r0=np.load(f,allow_pickle=False)



#%% Calculate Q
beta=2
Lambda=1
Q=(1/(1+np.exp(beta*(R-Lambda*r0))))

i=np.std(Q,axis=0)>.3
Q=Q.T[i].mean(0)

fig,ax=plt.subplots(3,1)
for k,a in enumerate(ax):
    a.plot(Q[len(Q)//3*k:len(Q)//3*(k+1)])