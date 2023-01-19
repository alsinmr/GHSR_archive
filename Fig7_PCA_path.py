#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 16:31:32 2022

@author: albertsmith
"""

import pyDR
from pyDR.PCA import PCA
import numpy as np
import os
import matplotlib.pyplot as plt
from pyDR.PCA.mouse_path import Path2Energy

# from pyDR.chimeraX.chimeraX_funs import set_chimera_path
# set_chimera_path() #Put your own path to ChimeraX here!!
#e.g. set_chimera_path('/Applications/ChimeraX-1.2.5.app/Contents/MacOS/ChimeraX')


#%% File locations
# mddir='/work/public/ghrelin-receptor'
mddir='/Volumes/My Book/GHSR'

topos=['WT-apo_run1_just_protein.pdb','WT-ghrelin_run1_just_protein.pdb']
states=['apo','bound']

#%% Re-create the PCA object from file

#String for selections in chimeraX
hlx_ind=np.concatenate((np.arange(42,69),np.arange(81,106),np.arange(116,143),
                        np.arange(162,183),np.arange(212,221),np.arange(225,234),
                        np.arange(257,287),np.arange(302,327),np.arange(327,340)))
hlx_sel='setattr :'+','.join(str(int(i)) for i in hlx_ind)+\
        ' res ss_type 1'
        
cmx_commands=['lighting soft','graphics silhouettes true',hlx_sel,'ribbon']

proj=pyDR.Project('Projects/PCA',create=True)
proj.chimera.saved_commands=cmx_commands

topo='WT-apo_run1_just_protein.pdb'
state='apo'


pca=PCA(pyDR.MolSelect(topo=os.path.join(mddir,topo),project=proj)).select_atoms('name N C CA')
with open(os.path.join('PCA_results',f'{state}_covar.data'),'rb') as f:
    pca._covar=np.load(f,allow_pickle=False)
with open(os.path.join('PCA_results',f'{state}_pcamp.data'),'rb') as f:
    pca._pcamp=np.load(f,allow_pickle=False)
with open(os.path.join('PCA_results',f'{state}_mean.data'),'rb') as f:
    pca._mean=np.load(f,allow_pickle=False)    
with open(os.path.join('PCA_results',f'{state}_Lambda.data'),'rb') as f:
    pca._lambda=np.load(f,allow_pickle=False) 
with open(os.path.join('PCA_results',f'{state}_PC.data'),'rb') as f:
    pca._PC=np.load(f,allow_pickle=False)        
p2e=Path2Energy(pca)

p2e.create_plots(cmap='nipy_spectral',maxbin=120,nmax=1,mode='points')
pts=np.array([[ -0.96774194,  79.48051948],  #These points define the path through the trajectory
       [ -1.4516129 ,  53.50649351],
       [ -3.38709677,  25.58441558],
       [-24.19354839,  12.5974026 ],
       [-45.48387097,  -0.38961039],
       [-65.32258065,  -9.48051948],
       [-50.80645161, -27.01298701],
       [-15.48387097, -14.67532468],
       [  0.        , -28.31168831],
       [ 17.90322581, -23.76623377],
       [ 29.03225806, -34.80519481],
       [ 41.61290323, -32.20779221],
       [ 59.51612903, -27.01298701],
       [ 86.61290323,  -8.18181818],
       [112.25806452,  17.79220779]]).T
p2e.load_points(pts)

p2e.plot_DelG(d=5)



indices=[1,6,8,10,13]  #These are the structures we'll actually show

p2e.chimera(pts=indices)

turns=['~show','turn x -90','view','turn y 90','turn z 15','turn y 45','turn z 5']
proj.chimera.command_line(turns)



cmap=plt.get_cmap('tab10')
color=[[str(int(c*100)) for c in cmap(k)] for k in [1,7,4,3,2]]
for c,k in zip(color,range(5)):
    proj.chimera.command_line(f'color #{k+1} '+','.join(c))
    
for m,view in enumerate([[],['turn x 90','zoom 1.5'],['turn x 180']]):
    proj.chimera.command_line(view)
    for n,k in enumerate([0,4,3,2]):
        proj.chimera.command_line(['~ribbon',f'ribbon #2|#{k+1}'])
        proj.chimera.savefig(f'apo_path_view{m+1}_pt{n}.png',options='transparentBackground true')
    




#%% Run detector analysis on first four principle components

from pyDR.misc.Averaging import avgDataObjs


n=7 #Number of detectors

lengths=[[372685,364134,328181],[364893,357953,342154]]

for q,(topo,state,pts0,length) in enumerate(zip(topos,states,pts,lengths)):
    if len(proj)<(q+1)*6:
        #Load the correct pca
        pca=PCA(pyDR.MolSelect(topo=os.path.join(mddir,topo),
                               project=proj)).select_atoms('name N C CA')
        with open(os.path.join('PCA_results',f'{state}_covar.data'),'rb') as f:
            pca._covar=np.load(f,allow_pickle=False)
        with open(os.path.join('PCA_results',f'{state}_pcamp.data'),'rb') as f:
            pca._pcamp=np.load(f,allow_pickle=False)
        with open(os.path.join('PCA_results',f'{state}_mean.data'),'rb') as f:
            pca._mean=np.load(f,allow_pickle=False)
        with open(os.path.join('PCA_results',f'{state}_Lambda.data'),'rb') as f:
            pca._lambda=np.load(f,allow_pickle=False) 
        with open(os.path.join('PCA_results',f'{state}_PC.data'),'rb') as f:
            pca._PC=np.load(f,allow_pickle=False)
        pca._t=np.linspace(0,3*355000*.1,3*355000)
            
        fit0=list()
        for k in range(3):
            data=pca.PCA2data(t0=sum(length[:k]),tf=sum(length[:k+1]))
            
            data.detect.r_no_opt(n)
            fit0.append(data.fit())
            print(k)
        proj.clear_memory()
        fit0=avgDataObjs(proj['no_opt'][-3:])
        fit0.detect.r_auto(n)
        fit=fit0.fit()
        fit=fit.opt2dist(rhoz_cleanup=True)
        
        proj.save()
        
#%% Plot the results
fig,ax=plt.subplots(10,2,squeeze=False)
clr=plt.get_cmap('tab10')
for q,(ax0,key) in enumerate(zip(ax.T,['apo','ghrelin'])):
    fit=proj['opt_fit']['.+'+key]
    ax0[0].set_title(key)
    for m,a in enumerate(ax0):
        if not(a.get_subplotspec().is_last_row()):a.set_xticklabels('')    
        for k in range(n):
            a.bar(k,fit.R[m,k],color=clr(k),edgecolor='black')


plt.show()

#%% Save the results to a text file
from misc_functions import save_detectors,save_PCA_path_energy
titles=['PCA components 0-9 for apo GHSR','PCA components 0-9 for ghrelin-bound GHSR']
save_detectors(7,[d for d in proj['opt_fit']],titles)
os.rename('source_data/Figure7.txt', 'source_data/Figure7_PCAdet.txt')

save_PCA_path_energy(7,p2e.get_path(), p2e.get_DelG(),title='Apo GHSR')
os.rename('source_data/Figure7.txt', 'source_data/Figure7_pathEnergy.txt')