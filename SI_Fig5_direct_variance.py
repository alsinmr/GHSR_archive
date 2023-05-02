#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 13:14:43 2023

@author: albertsmith
"""

import pyDR



projs=[pyDR.Project('Projects/backboneHN'),pyDR.Project('Projects/CACB')]

for proj in projs:
    #Fit the individual trajectories
    sub=proj['no_opt']-proj['AvOb']
    sub.detect.r_auto(7)
    sub.fit()
    sub=proj['proc']-proj['AvOb']
    sub.opt2dist(rhoz_cleanup=True)
    
    #Plot the results
    proj.close_fig('all')
    for m,key in enumerate(['apo','ghrelin']):
        proj.current_plot=m+1
        colors=[None,'black','grey']
        for k,color in enumerate(colors):
            if color is None:
                (proj['opt_fit']['.+'+key]-proj['AvOb'])[k].plot()
            else:
                po=(proj['opt_fit']['.+'+key]-proj['AvOb'])[k].plot(color=color)
        po.fig.set_size_inches([6.8,9.3])
        po.ax_sens.set_xlim([-11.5,-3.5])
        for a in po.ax:
            a.set_xlim([sub[0].label[0]-2,sub[0].label[-1]+2])
            a.set_ylim([0,1])
            
        proj.savefig('detector_var_'+key+'.pdf')