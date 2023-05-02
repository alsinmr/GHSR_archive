#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 12:49:04 2023

@author: albertsmith
"""



#%% Calculate Q
beta=5
Lambda=1.
Q=(1/(1+np.exp(beta*(R-Lambda*r0)))).mean(1)
ax=plt.figure().add_subplot(111)
ax.plot(Q)