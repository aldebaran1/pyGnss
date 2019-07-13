#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 14:50:28 2019

@author: smrak
"""

import xarray
import numpy as np
from datetime import datetime
import os
from pyGnss import pyGnss
from glob import glob
import georinex as grx
import matplotlib.pyplot as plt
from matplotlib import dates
from pymap3d import ecef2aer
from scipy.interpolate import CubicSpline

def sp3(f):
    with open(f, 'r') as F:
        # Read in the context
        sp3 = F.readlines()[22:-1]
        # Paramter initialization
        svlist = np.array([l[1:4] for l in sp3[1:33]])
        ecef = np.nan * np.ones((svlist.size,len(sp3[::33]), 3), dtype=np.float32)
        clock = np.nan * np.ones((svlist.size, len(sp3[::33])), dtype=np.float32)
        # Read in data
        # 1) Time
        navtimes = np.array([datetime.strptime(l.strip('\n')[3:-2], "%Y %m %d %H %M %S.%f") for l in sp3[::33]], dtype='datetime64[s]')
        for i, sv in enumerate(svlist):
            sp3 = sp3[1:]
            for j in range(navtimes.size):
                ecef[i,j,:] = np.array([x for x in sp3[::33][j][4:-1].lstrip().rstrip().split(' ') if x][:3], dtype=np.float32)
                clock[i,j] = np.float32(sp3[::33][0][47:60])
    F.close()

    nav = xarray.Dataset(coords={'time': navtimes, 'sv': svlist, 'xyz' : ['ecefx', 'ecefy', 'ecefz']})
    nav['ecef'] = (('sv', 'time', 'xyz'), ecef)
    nav['clock'] = (('sv', 'time'), clock)
    
    return nav

def gpsSatPositionSP3(fnav, dt, sv=None, rx_position=None, coords='xyz'):
    assert sv is not None
    # Read in data
    D = grx.load(fnav).sel(sv=sv)
    dt = dt.astype('datetime64[s]')
    navtimes = D.time.values.astype('datetime64[s]')
    CSx = CubicSpline(navtimes.astype(int), D.ecef.values[:,0])
    CSy = CubicSpline(navtimes.astype(int), D.ecef.values[:,1])
    CSz = CubicSpline(navtimes.astype(int), D.ecef.values[:,2])
    ecefxi = CSx(obstimes.astype(int))
    ecefyi = CSy(obstimes.astype(int))
    ecefzi = CSz(obstimes.astype(int))
    
    if coords == 'xyz':
        return np.array([ecefxi, ecefyi, ecefzi])
    else:
        AER = ecef2aer(x=ecefxi, y=ecefyi, z=ecefzi,
           lon0=rx_position[1], lat0=rx_position[0], h0=rx_position[2])
        return AER



fnc = 'morm2440.17d.nc'
fnav = 'brdc2440.17n'
DOBS = pyGnss.dataFromNC(fnc,fnav,sv='G07',tlim=None,el_mask=30, satpos=True)
obstimes = DOBS.time.values.astype('datetime64[s]')
RX = DOBS.position_geodetic

fsp3 = glob(os.getcwd() + '/*.sp3')[0]
D = grx.load(fsp3).sel(sv='G07')
navtimes = D.time.values.astype('datetime64[s]')
AER = ecef2aer(x=D.ecef.values[:,0], y=D.ecef.values[:,1], z=D.ecef.values[:,2],
           lon0=RX[1], lat0=RX[0], h0=RX[2])

ecef = pyGnss.gpsSatPositionSP3(fsp3, dt=obstimes, sv='G07', rx_position=RX, coords='xyz')
AERi = pyGnss.gpsSatPositionSP3(fsp3, dt=obstimes, sv='G07', rx_position=RX, coords='aer')

fig = plt.figure(figsize=[8,6])
ax0 = fig.add_subplot(211)
ax0.plot(D.time.values, D.ecef.values[:,0], '.b', label='ecefx')
ax0.plot(obstimes, ecef[0,:], 'b')
ax0.plot(D.time.values, D.ecef.values[:,1], '.r', label='ecefy')
ax0.plot(obstimes, ecef[1,:], 'r')
ax0.plot(D.time.values, D.ecef.values[:,2], '.k', label='ecefz')
ax0.plot(obstimes, ecef[2,:], 'k')
ax0.legend()
ax1 = fig.add_subplot(212, sharex=ax0)
ax1.plot(D.time.values, AER[0], '.b', label='azimuth')
ax1.plot(obstimes, AERi[0], 'b')
ax11 = ax1.twinx()
ax11.plot(D.time.values, AER[1], '.r', label='elevation')
ax11.plot(obstimes, AERi[1], 'r')
ax1.legend()
ax11.legend()
myFmt = dates.DateFormatter('%H:%M')
ax1.xaxis.set_major_formatter(myFmt)