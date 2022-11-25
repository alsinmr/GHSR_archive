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



import pyDR


#%% Load the project and fit/optimize the detectors
proj=pyDR.Project('Projects/backboneHN_half')

#We'll check if the operations have already been run/saved before executing them

if not(len(proj['opt_fit'])):
    proj['no_opt'].detect.r_auto(7)
    proj['no_opt'].fit()
    proj['proc'].opt2dist(rhoz_cleanup=True)
    proj.save()

#%% Plot the results

plot=0
for state in ['apo','ghrelin']:
    for k in range(2):
        plot+=1
        proj.current_plot=plot
        if state=='ghrelin' and k==1:
            proj['opt_fit']['.+'+state][[1,5]].plot()
        else:
            proj['opt_fit']['.+'+state][k::2].plot()
        proj.plot_obj.show_tc()
        title=f'{state.capitalize()}, '+('second' if k else 'first')+' half'
        proj.plot_obj.ax_sens.set_title(title)
        proj.plot_obj.fig.set_size_inches([6.7,10.3])
        proj.savefig(title.replace(',','').replace(' ','_'))
