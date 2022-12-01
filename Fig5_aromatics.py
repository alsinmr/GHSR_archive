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
import os

proj=pyDR.Project('Projects/aromatic_iRED')

#%% FInish processing
if not(len(proj['opt_fit'])):
    proj.detect.r_auto(7)
    proj.fit()['proc'].opt2dist(rhoz_cleanup=True)
    proj.save()
    
proj['opt_fit'].modes2bonds()

avgDataObjs(proj['iREDbond']['.+apo'])
avgDataObjs(proj['iREDbond']['.+ghrelin'])

sub=proj['iREDbond']['AvOb_rk1']



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
ax.set_xlabel('Detector')
ax.set_ylabel(r'$\rho_n^{(\theta,S)}$')
        
#%% Plot correlation coefficient for all detectors

fig=plt.figure()
fig.clear()
fig.set_size_inches([6.5,12])
ax0=fig.add_subplot(nd+1,1,1)
ax0.set_title('Correlation Coefficient')
ax=[fig.add_subplot(nd+1,1,k+2) for k in range(nd)]
sub[0].sens.plot_rhoz(index=range(nd),ax=ax0)

# for a in ax[1:]:a.sharex(ax[0])

AA=pyDR.tools.AA

for k,a in enumerate(ax):
    for m,(i0,i,s,h) in enumerate(zip(index0,index,sub,hatch)):
        x=np.arange(sum(i))+w/2*(-1)**(m+1)
        a.bar(x,np.abs(s.CCnorm[k][i0][i]),w,color=cmap(k),hatch=h)
        a.set_ylim([0,0.6])
        a.set_ylabel(r'$\rho_n^{CC}$')
        a.set_yticks(np.linspace(0,.6,4))
        
    a.set_xticks(np.arange(sum(i)))

    if a.is_last_row():
        a.set_xticklabels([str(sel[0].resid)+AA[sel[0].resname].symbol for sel in s.select.sel1[i]],rotation=90)
    else:
        a.set_xticklabels([])

#%% Plot Chemical Shift prediction from SHIFTX2
def load_file(filename='SHIFTX2/apo_run1_CA.txt'):
    out={}
    with open(filename,'r') as f:
        for line in f:
            resid,*shifts=line.strip().split('\t')
            out[int(resid)]=np.array(shifts).astype(float)
    return out

apoCA=np.zeros([len(resids),3,710])
for k in range(3):
    out=load_file(f'SHIFTX2/apo_run{k+1}_CA.txt')
    for m,resid in enumerate(resids):
        apoCA[m,k]=out[resid]
        
boundCA=np.zeros([len(resids),3,710])
for k in range(3):
    out=load_file(f'SHIFTX2/ghrelin_run{k+1}_CA.txt')
    for m,resid in enumerate(resids):
        boundCA[m,k]=out[resid]
        

w=0.4
ax=plt.figure().add_subplot(111)
shift=apoCA.reshape([len(resids),3*710])
ax.bar(np.arange(len(resids))-w/2,shift.mean(1),color='grey',edgecolor='black',width=w)
# ax.errorbar(np.arange(len(resids))-w/2,shift.mean(1),shift.std(1),color='black',capsize=4,linestyle='')        
ax.errorbar(np.arange(len(resids))-w/2,shift.mean(1),apoCA.mean(-1).std(1),color='black',capsize=4,linestyle='')        
shift=boundCA.reshape([len(resids),3*710])
ax.bar(np.arange(len(resids))+w/2,shift.mean(1),color='darkgrey',edgecolor='black',width=w)
ax.errorbar(np.arange(len(resids))+w/2,shift.mean(1),boundCA.mean(-1).std(1),color='black',capsize=4,linestyle='')
ax.set_xticks(range(len(resids)))
ax.set_xticklabels(resids,rotation=90)
ax.set_ylabel(r'$\delta^{C\alpha}$ / ppm')
ax.set_ylim([56.5,62.5])


#%% Make chimera plots
res1=[*resids,276]
index=[np.isin(s.label,[str(r) for r in res1]) for s in sub]  #These are the residues we'll correlate with
sel=[copy(s.select) for s in sub]
for s,i in zip(sel,index):
    for k in ['repr_sel','sel1','sel2']:
        setattr(s,k,[x for x in getattr(s,k)[i]])

cmx_cmds=['~show','ribbon','show :'+','.join([str(r) for r in res1]),
          'view','turn x -90','turn z -5','turn x 5','align #2/B@CA toAtoms #1@CA',
          'color light steel blue target r','move x 50 models #2','zoom 1.2',
          'lighting full','graphics silhouettes True','color :276 30,30,30','view','zoom 1.15']

for k in range(6):
    sub.chimera.close()
    for sel0,s,i,i0 in zip(sel,sub,index,index0):
        x=s.CCnorm[k][i0][i]
        x*=1.75
        x[np.argwhere(np.argwhere(i)[:,0]==i0)[0,0]]=1.2
        sel0.chimera(x=x.T,color=cmap(k))
    sub.chimera.command_line(cmx_cmds)
    sub.chimera.savefig('toggle_switchCC_rho{0}.png'.format(k),options='transparentBackground True')
