#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:56:23 2022

@author: albertsmith
"""

import numpy as np


#%% Defines where helices are in GHSR
def load_helices():
    name=list()
    resids=list()
    with open('helix_assign.txt','r') as f:
        for line in f:
            name.append(line.strip().split(':')[0])
            ar=np.arange(\
                int(line.split('resid')[1].split('to')[0].strip()),
                int(line.split('to')[1].strip())+1)
            if len(resids):
                ar=ar[np.logical_not(np.isin(ar,np.concatenate(resids)))]
            resids.append(ar)
    return name,resids


def helix_only(select=None):
    names,hlx=load_helices()
    i=np.arange(32,340)
    i0=np.concatenate(np.array(hlx)[np.array(['L' in name or 'NT' in name for name in names],dtype=bool)])
    i0=np.logical_not(np.isin(i,i0))
    i=i[i0]
    if select is None:return i
    return np.array([s.resids[0] in i for s in select.sel1],dtype=bool)


import os

def save_detectors(fignum:int,data,titles,CCindex:list=None):
    """
    Saves detector or cross-correlation information to a text file.
    
    Parameters
    ----------
    fignum : 
        Figure number to which data corresponds (determines file name).
    data : list
        List of data objects that are being saved
    titles : list
        List of titles corresponding to each data object
    CCindex : list
        If storing CC data, these are the labels of the residues to correlate to

    Returns
    -------
    None.

    """
    
    #Function for writing out series of detector responses
    def write_responses(f,R,label):
        f.write(' '*10+''.join([f'{i:<9}' for i in range(R.shape[1])])+'\n')
        # f.write(' '*8+('R'+' '*7)*R.shape[1]+'\n')
        for R0,label0 in zip(R,label):
            label0=str(label0)
            f.write(f'{label0[:9]:>9s} '+' '.join([f'{np.abs(R00):8f}' for R00 in R0])+'\n')
    
    # File location
    folder='source_data/'
    filename=os.path.join(folder,f'Figure{fignum}.txt')
    
    
    if not(isinstance(data,list)):data=[data]   #Make sure a list
    if not(isinstance(titles,list)):titles=[titles]

    with open(filename,'w') as f:
        f.write('#'*60+'\n')
        f.write('-'*26+f'Figure {fignum}'+'-'*26+'\n')
        f.write('#'*60+'\n\n')

        for k,(data0,title) in enumerate(zip(data,titles)):
            f.write(f'\nDATA SET #{k+1}: {title}')
            if CCindex is None:
                f.write('\nDetector responses\n')
                write_responses(f,data0.R,data0.label)
            else:
                if not(isinstance(CCindex,list)):CCindex=[CCindex]
                for CC0 in CCindex:
                    f.write(f'\nDetector responses for {CC0}\n')
                    i=data0.label.tolist().index(CC0)
                    write_responses(f,data0.R[i:i+1],data0.label[i:i+1])
                    f.write(f'\nCross-correlation coefficients to {CC0}\n')
                    write_responses(f,data0.CCnorm[:,i].T,data0.label)
        
     
def save_PCA_hist(fignum:int,pca:list,titles:list,nbins=100,maxbins=120):
    # File location
    folder='source_data/'
    filename=os.path.join(folder,f'Figure{fignum}.txt')

    def write_hist(f,pca,nbins,maxbins):
        bins=np.linspace(-maxbins,maxbins,nbins)
        for n0 in range(3):
            f.write(f'\nHistogram between components {n0} and {n0+1}\n')
            f.write(' '*6+' '.join([f'{b:>6.1f}' for b in bins])+'\n')
            hist=np.histogram2d(pca.PCamp[n0],pca.PCamp[n0+1],bins=bins)[0].astype(int)
            for b,h in zip(bins,hist):
                f.write(f'{b:>6.1f}'+' '.join([f'{h0:6}' for h0 in h])+'\n')
        
    
    if not(isinstance(titles,list)):titles=[titles]
    if not(isinstance(pca,list)):pca=[pca]
    
    with open(filename,'w') as f:
        f.write('#'*60+'\n')
        f.write('-'*26+f'Figure {fignum}'+'-'*26+'\n')
        f.write('#'*60+'\n\n')
        
        for k,(pca0,title) in enumerate(zip(pca,titles)):
            f.write(f'\nDATA SET #{k+1}: Histogram for {title}\n')
            write_hist(f,pca0,nbins=nbins,maxbins=maxbins)
            
def save_PCA_path_energy(fignum:int,path,DelG,title=None):
    # File location
    folder='source_data/'
    filename=os.path.join(folder,f'Figure{fignum}.txt')
    
    with open(filename,'w') as f:
        if title is not None:f.write(title+'\n')
        f.write('Delta G along path\n')
        f.write('       PC0       PC1      DelG\n')
        for (PC0,PC1),delG in zip(path.T,DelG):
            f.write(f'{PC0:>10.3f}{PC1:>10.3f}{delG:>10.3f}\n')
    