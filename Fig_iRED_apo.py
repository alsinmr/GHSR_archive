#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 21:25:28 2022

@author: albertsmith
"""

import pyDR
import numpy as np
import matplotlib.pyplot as plt
from pyDR.misc.Averaging import avgDataObjs

proj=pyDR.Project('Projects/apo_iRED')

if not(len(proj['opt_fit'])):
    proj.detect.r_auto(7)
    proj.fit()['proc'].opt2dist(rhoz_cleanup=True)
    proj.save()


proj['opt_fit'].modes2bonds()
avgDataObjs(proj['iREDbond'])


data=proj['iREDbond'][-1]