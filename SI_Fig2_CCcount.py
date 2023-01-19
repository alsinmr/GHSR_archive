#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 14:46:46 2023

@author: albertsmith
"""

import pyDR
import numpy as np
import matplotlib.pyplot as plt
from pyDR.misc.Averaging import avgDataObjs

"If first time running pyDR with chimeraX, uncomment the lines below, and provide chimeraX path"
# from pyDR.chimeraX.chimeraX_funs import set_chimera_path
# set_chimera_path() #Put your own path to ChimeraX here!!
"e.g. set_chimera_path('/Applications/ChimeraX-1.2.5.app/Contents/MacOS/ChimeraX')"

proj=pyDR.Project('Projects/ghrelin_iRED')


proj['opt_fit'].modes2bonds()
avgDataObjs(proj['iREDbond'])

data=proj['iREDbond'][-1]


#%% Interactions with correlation over 0.5
#(SI figure 1)

hlx=[np.arange(41,71),np.arange(78,106),np.arange(114,146),np.arange(158,182),
     np.arange(209,238),np.arange(253,290),np.arange(296,327),np.arange(327,340)]
ter=[np.arange(33,41)]
icl=[np.arange(71,78),np.arange(146,158),np.arange(238,253),]
ecl=[np.arange(106,114),np.arange(182,209),np.arange(290,296)]

groups=[ter[0],hlx[0],icl[0],hlx[1],ecl[0],hlx[2],icl[1],hlx[3],ecl[1],hlx[4],
        icl[2],hlx[5],ecl[2],hlx[6],hlx[7]]
names=['N-term','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','H8']


resids=[int(lbl.split('_')[0]) for lbl in data.label]

cutoff=0.5

count1to7=list()
count9to15=list()
for g in groups:
    index=np.isin(resids,g)
    count1to7.append((np.abs(data.CCnorm[5][:5][:,index])>cutoff).sum())
    count9to15.append((np.abs(data.CCnorm[5][5:13][:,index])>cutoff).sum())

fig=plt.figure()
ax=[fig.add_subplot(2,1,k+1) for k in range(2)]

for a,count in zip(ax,[count1to7,count9to15]):
    a.bar(range(len(count)),count,edgecolor='black',width=.6)
    a.set_ylim([0,max(*count1to7,*count9to15)+1])
    a.set_xticks(range(len(count)))
    a.set_xticklabels(names,rotation=90)
    a.set_ylabel('Count')

ax[0].legend((r'ghrelin$_{2-6}$',))
ax[1].legend((r'ghrelin$_{8-15}$',))
ax[0].set_xticklabels([])
fig.tight_layout()    
        
plt.show()