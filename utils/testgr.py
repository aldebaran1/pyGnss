#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 17:08:05 2019

@author: smrak
"""
from pyGnss import pyGnss
from datetime import datetime
import matplotlib.pyplot as plt
sv = 'G12'
fnav = '/media/smrak/gnss/nav/brdc0460.13n'
fnc = '/media/smrak/gnss/obs/2013/046/ztl40460.13d.nc'
el_mask = 30

tlim = [datetime(2013, 2, 15, 11, 0, 0),
        datetime(2013, 2, 15, 13, 0, 0)]
#tlim = None

D = pyGnss.dataFromNC(fnc,fnav,sv=sv,tlim=tlim,el_mask=el_mask, 
                      satpos=True, ipp=True, ipp_alt = 110)

plt.plot(D.time.values, D['S1'].values)