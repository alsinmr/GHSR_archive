#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 14:48:37 2022

@author: albertsmith
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:15:16 2022

@author: albertsmith
"""


import sys
sys.path.append('/work/home/alsi/GitHub')
import pyDR
from pyDR.misc.Averaging import avgDataObjs


#%% Load the project and fit/optimize the detectors
proj=pyDR.Project('Projects/backboneHN_half')

#We'll check if the operations have already been run/saved before executing them

if not(len(proj['opt_fit'])):
    proj['no_opt'].detect.r_auto(7)
    proj['no_opt'].fit()['proc'].opt2dist(rhoz_cleanup=True)
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
      'turn x -90','turn y 15 models #1','turn y 15 models #2',
      'lighting full','graphics silhouettes true','view','zoom 1.15']
for k in range(proj[-1].R.shape[1]):
    proj.chimera.close()
    proj['opt_fit'].chimera(scaling=1.5,rho_index=k)
    proj.chimera.command_line(cmds)
    proj.chimera.savefig('rho{0}'.format(k),'transparentBackground True')