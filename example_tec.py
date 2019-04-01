#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 11:55:08 2019

@author: smrak
"""

import numpy as np
import georinex as gr
from pyGnss import pyGnss
from pyGnss import gnssUtils as gu
import matplotlib.pyplot as plt

def _chi2(y):
    cost = np.nan * np.ones(y.shape[1])
    for it in range(y.shape[1]):
        idx = np.where(np.isfinite(y[:,it]))[0]
        d = []
        for i in range(idx.shape[0]):
            for j in np.arange(i+1, idx.shape[0]):
                d.append(abs(y[idx[i], it] - y[idx[j], it]))
        cost[it] = np.nanmean(d)
        
    return cost

# Files
NCF = '/media/smrak/gnss/obs/2017/147/mojc1470.17d.nc'
NF  = '/media/smrak/gnss/nav/brdc1470.17n'
SBF = '/media/smrak/gnss/jplg/jplg1470.17i'

SB = True
el_mask = 30
tlim = None
H = 350

svlist = gr.load(NCF).sv.values
navdata = gr.load(NF)
navdatatime = navdata.time.values
D = gr.load(NCF)
time = D.time.values
satbias = pyGnss.getSatBias(SBF)

fig1 = plt.figure(figsize = [8,6])
fig2 = plt.figure(figsize = [8,6])
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)

stec = np.nan * np.ones((svlist.shape[0], D.time.values.shape[0]))
F = np.nan * np.ones((svlist.shape[0], D.time.values.shape[0]))
for i, sv in enumerate(svlist):
    D = pyGnss.dataFromNC(NCF, NF, sv=sv, tlim=tlim, el_mask=el_mask, satpos=True, 
                          ipp=True, ipp_alt = H)
    idel = D['idel'].values
    sb = satbias[sv]
    
    dt = D.time.values
    tsps = np.diff(dt.astype('datetime64[s]'))[0].astype(int)
    el = D.el.values
    C1 = D['C1'].values
    C1[~idel] = np.nan
    C2 = D['P2'].values
    C2[~idel] = np.nan
    L1 = D['L1'].values
    L1[~idel] = np.nan
    L2 = D['L2'].values
    L2[~idel] = np.nan
    
#    obs = {'C1': C1, 'P2': C2, 'L1': L1, 'L2': L2}
    
#    tec, tecd, F = pyGnss.processTEC(obs, sv, frequency = 2, Ts=tsps, 
#                                     H=H, elevation = el, sat_bias=sb)
    
    # Intervals
    ixin, intervals = pyGnss.getIntervals(L1,L2,C1,C2, maxgap=1)
    tec = np.nan * stec[0]
    for r in intervals:
        tec[r[0]:r[-1]] = pyGnss.slantTEC(C1[r[0]:r[-1]], C2[r[0]:r[-1]], 
                                          L1[r[0]:r[-1]], L2[r[0]:r[-1]])
    stec[i, :] = tec + sb
    F[i, :] = pyGnss.getMappingFunction(el, H)

vtec = stec * F

rxb = -5
eps = 1e6
#cost_column0 = np.nanmean(vtec[:,0])
#cost0 = np.nanmean(cost_column0)
#cost = np.nan * np.ones(stec.shape[1])
#for it in range(stec.shape[1]):
#    idx = np.where(np.isfinite(vtec[:,it]))[0]
#    d = []
#    for i in range(idx.shape[0]):
#        for j in np.arange(i+1, idx.shape[0]):
#            d.append(abs(vtec[idx[i], it] - vtec[idx[j], it]))
#    cost[it] = np.nanmean(d)
CX = []
cost1 = _chi2(vtec)
CX.append(np.nansum(cost1))
while eps > 3:
    stec1 = np.add(np.copy(stec), rxb)
    vtec1 = np.multiply(stec1, F)
    cost2 = _chi2(vtec1)
    CX.append(np.nansum(cost2))
    if CX[-1] < CX[-2]:
        continue
    elif cost2 > cost1:
        
#    cost_column2 = np.divide(np.nansum(vtec1, axis = 0), np.nanmean(vtec1, axis = 0))
#    cost2 = np.nanmean(cost_column2)
    
    break
    
fig1 = plt.figure(figsize = [8,6])
ax1 = fig1.add_subplot(111)
for y in vtec:
    ax1.plot(time, y)
    
fig1 = plt.figure(figsize = [8,6])
ax1 = fig1.add_subplot(111)
for y in vtec1:
    ax1.plot(time, y)
    
plt.figure(figsize = [6,3])
plt.plot(time, cost1, 'b')
plt.plot(time, cost2, 'r')


#    tecd = gu.getPlainResidual(stec, Ts=tsps)
#    F = pyGnss.getMappingFunction(el, H)
    
#    while eps > 10:
#        stec + rxb
#    ax1.plot(dt, stec)
#    ax2.plot(dt, stec * F)