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
from misc_functions import color_code

#%% Function for color-coding protein regions



#%% Load the project and optimize/fit the detectors
proj=pyDR.Project('Projects/backboneHN')
proj['no_opt'].detect.r_auto(7)
proj['no_opt'].fit()
proj['proc'].opt2dist(rhoz_cleanup=True)


avgDataObjs(proj['.+WT-apo']['opt_fit'])
avgDataObjs(proj['.+WT-ghrelin']['opt_fit'])

proj['opt_fit']['.+AvOb'].plot()
for a in proj.plot_obj.ax:a.set_ylim([0,1])
proj.plot_obj.ax[0].set_xlim([30,340])
proj.plot_obj.show_tc()
proj.plot_obj.ax_sens.set_xlim([-11.5,-3.5])
proj.plot_obj.fig.set_size_inches([6.53, 9.01])

proj.savefig('apo_v_bound.pdf')