#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:15:16 2022

@author: albertsmith
"""


import pyDR
from pyDR.misc.Averaging import avgDataObjs
import matplotlib.pyplot as plt

"If first time running pyDR with chimeraX, uncomment the lines below, and provide chimeraX path"
# from pyDR.chimeraX.chimeraX_funs import set_chimera_path
# set_chimera_path() #Put your own path to ChimeraX here!!
"e.g. set_chimera_path('/Applications/ChimeraX-1.2.5.app/Contents/MacOS/ChimeraX')"


#%% Load the project and fit/optimize the detectors
proj=pyDR.Project('Projects/backboneHN')

#We'll check if the operations have already been run/saved before executing them

#Average the 3x2 trajectories
if not(len(proj['.+AvOb'])):
    avgDataObjs(proj['.+apo'])
    avgDataObjs(proj['.+ghrelin'])

if not(len(proj['opt_fit'])):  #Convert no-opt detectors to 7 optimized detectors
    proj['no_opt']['.+AvOb'].detect.r_auto(7)
    proj['no_opt']['.+AvOb'].fit()
    proj['proc']['.+AvOb'].opt2dist(rhoz_cleanup=True)
    proj.save()

#%% Plot the results
proj['opt_fit']['.+AvOb'].plot()

for a in proj.plot_obj.ax:a.set_ylim([0,1])
proj.plot_obj.ax[0].set_xlim([30,340])
proj.plot_obj.show_tc()
proj.plot_obj.ax_sens.set_xlim([-11.5,-3.5])
proj.plot_obj.fig.set_size_inches([6.53, 9.01])

proj.savefig('apo_v_bound.pdf')

proj.chimera.saved_commands=[]
cmds=['~show ~@N,C,CA,H','~show #2/A','show #2/A@N,C,CA','color #2/A slate grey',
      'turn x -90','turn y 15 models #1','turn y 15 models #2',''
      'lighting full','graphics silhouettes true','view','zoom 1.15',
      f'move x -{proj[-1].select.box[0]} coordinateSystem #2 atoms #2/A']
for k in range(proj[-1].R.shape[1]):
    proj.chimera.close()
    proj['opt_fit'].chimera(scaling=1.5,rho_index=k)
    proj.chimera.command_line(cmds)
    proj.chimera.savefig('rho{0}'.format(k),'transparentBackground True')
    
plt.show()

#%% Save results to a text file (attached to paper, also in github in source_data folder)
from misc_functions import save_detectors

data=[d for d in proj['opt_fit']]  #Needs to be a list
titles=['H–N motion for apo GHSR','H–N motion for ghrelin-bound GHSR']
save_detectors(fignum=1,data=data,titles=titles)