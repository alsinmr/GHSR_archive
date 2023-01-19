#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 15:00:31 2022

@author: albertsmith
"""



import pyDR
from pyDR.PCA import PCA
import numpy as np
import os
import matplotlib.pyplot as plt

"If first time running pyDR with chimeraX, uncomment the lines below, and provide chimeraX path"
# from pyDR.chimeraX.chimeraX_funs import set_chimera_path
# set_chimera_path() #Put your own path to ChimeraX here!!
"e.g. set_chimera_path('/Applications/ChimeraX-1.2.5.app/Contents/MacOS/ChimeraX')"

#%% File locations
# mddir='/work/public/ghrelin-receptor'
mddir='/Volumes/My Book/GHSR'

topos=['WT-apo_run1_just_protein.pdb','WT-ghrelin_run1_just_protein.pdb']
states=['apo','bound']



        

#%% Setup for using chimera for PCA display
proj=pyDR.Project('Projects/PCA',create=True)  #We'll just use project for saving some figures
chimera=proj.chimera       #Get the chimeraX manager from the project
hlx_ind=np.concatenate((np.arange(42,69),np.arange(81,106),np.arange(116,143),
                        np.arange(162,183),np.arange(212,221),np.arange(225,234),
                        np.arange(257,287),np.arange(302,327),np.arange(327,340)))  #These residues will be displayed as helices
hlx_sel='setattr :'+','.join(str(int(i)) for i in hlx_ind)+\
        ' res ss_type 1'
cmx_commands=['lighting soft','graphics silhouettes true',hlx_sel,'ribbon']  #Commands to pass to ChimeraX


#These pts correspond to 4-5 structures shown in Figure 6
#For each point, we have to define the location in the first 4 principal components
pts=[np.array([[-50.14460887, -26.77804288,  38.40141356,  -7.17965505],
       [ -3.81920641,  54.73120956,  -3.23280891, -11.08068449],
       [ 86.48600853,  -9.1861179 ,  -2.06001391,  -4.57896876],
       [ 30.77824607, -35.57400538, -28.44790139,  11.025149  ],
       [  0.28557609, -29.71003038,   8.49514109,  -7.17965505]]).T,
     np.array([[-31.01165213,  -2.00941431,  30.32641408,  -1.71542553],
            [ 64.32903773,  -2.67613242,   3.65768964,  -2.43351064],
            [-11.67682692,  51.32803457, -26.34462536,   5.46542553],
            [-28.34477969, -27.34470252, -35.67867891,  -8.8962766 ]]).T]




pca0=list()
for n,(topo,state,pts0) in enumerate(zip(topos,states,pts)):
    """
    Note that the pyDR PCA module is not fully set up for export, so this is
    a little hack-y. We have saved parts of the analysis on a server and
    reload them into the hidden variables in pca to restore functionality.
    """

    chimera.close()
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
    pca._t=np.linspace(0,3*355000*.1,3*355000)
    pca0.append(pca)
        
            
    # w,v=eigsh(pca.CoVar,k=10,which='LM')    
    # i=np.argsort(w)[::-1]
    # pca._lambda,pca._PC=w[i],v[:,i]

    "Plot the 0th-3rd components of the PCA"
    #This creates an interactive plot, from which Figure 6 was generated      
    pca.project.chimera.saved_commands=cmx_commands
    fig=pca.hist2struct(nmax=3,maxbin=120,cmap='nipy_spectral',nbins=100,ref_struct=False)[0].figure
    #From here we could interactively select the points to make the plots (just click on the plots)
    #However, in this case, we will just load the saved points from above
    pca.load_points(pts0) #This will mark the points in the PCA histograms and open the structures in ChimeraX
    
    fig.set_size_inches([9,8])
    fig.tight_layout()
    
    if state==states[0]:
        chimera.command_line('color #1 grey')
        plt.show() #Script will pause here until the figure is closed!
    else:
        chimera.command_line(['color #1/A blue','color #2/A yellow','color #3/A darkred','color #4/A darkgreen'])
        plt.show()
    #Use chimera.savefig(filename) to save views in chimera to the project folder
    


            
#%% Save histogram date to a text file
from misc_functions import save_PCA_hist
titles=['Histogram for PCA components for Apo GHSR, Components(0,1),(1,2),(2,3)',
        'Histogram for PCA components for Ghrelin-bound GHSR, Components(0,1),(1,2),(2,3)']
save_PCA_hist(6, pca0, titles)            
        
    
    
    

    
        

