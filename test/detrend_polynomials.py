# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 16:19:39 2018

@author: smrak
"""
from datetime import datetime, timedelta
from pandas import Timestamp
import numpy as np
from pyGnss import pyGnss
from pyGnss import gnssUtils as gu
import georinex as gr
import matplotlib.pyplot as plt

def mahali(folder='E:\\mahali\\',day=280,rx=9):
    # Init
    date = datetime.strptime(str(2015)+str(day), '%Y%j')
    directory = date.strftime('%Y-%m-%d')
    rxmah = {'2': 'mah2'+str(day)+'0.15o.nc',
             '3': 'mah3'+str(day)+'0.15o.nc',
             '4': 'mah4'+str(day)+'0.15o.nc',
             '5': 'mah5'+str(day)+'0.15o.nc',
             '6': 'mah6'+str(day)+'0.15o.nc',
             '7': 'mah7'+str(day)+'0.15o.nc',
             '8': 'mah8'+str(day)+'0.15o.nc',
             '9': 'mah9'+str(day)+'0.15o.nc',
             '13': 'ma13'+str(day)+'0.15o.nc',}
    folder = folder + directory + '\\'
    if isinstance(rx, int):
        nc = rxmah[str(rx)]
        fnc = folder + nc
    elif isinstance(rx, list):
        fnc = []
        for r in rx:
            nc = rxmah[str(rx)]
            fnc.append(folder + nc)
    nav = 'brdc'+str(day)+'0.15n'
    fnav = folder + nav
    
    return fnc, fnav

def detrendMetric(x, porders=[1,15]):
    value = []
    totalsum = []
    for po in porders:
        ptemp = pyGnss.phaseDetrend(x,po)
        value.append(np.sum(ptemp))
        totalsum.append(np.sum(abs(ptemp)))
    return np.array(value), np.array(totalsum)

def dataFromNC(fnc,fnav,el_mask=30,tlim=None):
    leap_seconds = gu.getLeapSeconds(fnav)
    D = gr.load(fnc, useindicators=True).sel(sv=sv)
    if tlim is not None:
        if len(tlim) == 2:
            D = D.where(np.logical_and(D.time >= np.datetime64(tlim[0]), D.time <= np.datetime64(tlim[1])), drop=True)
    obstimes64 = D.time.values
    dt = np.array([Timestamp(t).to_pydatetime() for t in obstimes64]) - \
               timedelta(seconds = leap_seconds)
    rx_xyz = D.position
    aer = pyGnss.gpsSatPosition(fnav,dt,sv=sv, rx_position=rx_xyz, coords='aer')
    idel = (aer[1] >= el_mask)
    dt = dt[idel]
    return dt, D, idel

def getDoppler(x):
    D = np.diff(x) # Diff 1st order
    D = np.hstack((np.nan, D))
    return D

folder = 'C:\\Users\\smrak\\Google Drive\\BU\\projects\\Themahali\\data\\rinex\\'
sv = 'G02'
rx = 7
day = [280,288]
year = 2015
plt.figure(figsize=(8,5))
plt.title('SV: '+sv)
for day in np.arange(day[0],day[1]+1):

    tlim = [datetime.strptime("{} {} {}:{}".format(year,day,6,5), "%Y %j %H:%M"), 
            datetime.strptime("{} {} {}:{}".format(year,day,12,5), "%Y %j %H:%M")]
    fnc, fnav = mahali(folder = folder, rx = rx, day = day)
    
    dt, D, idel = dataFromNC(fnc,fnav,tlim=tlim)
    
    # L1
    L1 = D['C1'].values[idel]
    prange = [2,15]
    porders = np.arange(prange[0], prange[1]+1)
    mean, residualsum = detrendMetric(L1, porders)

#plt.figure(figsize=(8,5))
    plt.semilogy(porders, residualsum, '.', markersize=5)
    plt.semilogy(porders, residualsum, '--', label=str(day))
    plt.xlabel('Polynomial order')
    plt.ylabel('Residual sum')
plt.legend()