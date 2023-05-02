#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 15:59:41 2023

@author: albertsmith
"""

import pyDR

proj=pyDR.Project('Projects/backboneHN')


#%% Estimate R2 from the MD trajectories
nmr=pyDR.Sens.NMR(Type='R1p',v1=0,vr=12.777,v0=600,Nuc='13CA') #R2 sensitivity for Calpha

sub=proj['no_opt']['AvOb']      #Starting point for estimating R2
n=10
sub.detect.r_target(nmr,n=n)   #Set detector sensitivity to approximate R2 sensitivity

sub.fit()   #Solve the R2 values

ax=sub[0].detect.plot_rhoz(index=[0])[0].axes    #Plot sensitivities
nmr.plot_Rz(ax=ax,color='black',linestyle=':')


i_apo=[lbl in [186,258,280] for lbl in sub['.+apo'].label]
i_bound=[lbl in [186,258,280] for lbl in sub['.+ghrelin'].label]

DelR2=proj[f'p{n}.+apo'].R[i_apo][:,0]-proj[f'p{n}.+ghrelin'].R[i_bound,0]

Del_ppm=DelR2/150/np.pi