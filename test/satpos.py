# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 11:28:03 2018

@author: smrak
"""

import numpy as np
import georinex as gr
import os
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from pyGnss import pyGnss
import xarray
from pandas import Timestamp
from pymap3d import ecef2aer, ecef2geodetic

def gpsSatPosition(fnav, dt, sv=None):
    navdata = gr.load(fnav).sel(sv=sv)
    timesarray = np.asarray(dt,dtype='datetime64[ns]') #[datetime64 [ns]]
    # Manipulate with times, epochs and crap like this
    navtimes = navdata.time.values # [datetime64 [ns]]
    idnan = np.isfinite(navdata['Toe'].values)
    navtimes = navtimes[idnan]
    bestephind = []
    for t in timesarray:
        idt = abs(navtimes - t).argmin() if t>navtimes[abs(navtimes - t).argmin()] else abs(navtimes - t).argmin()-1
        bestephind.append(idt)
#    bestephind = np.array([np.argmin(abs(navtimes-t)) for t in timesarray])
    gpstime = np.array([pyGnss.getGpsTime(t) for t in dt])
    t = gpstime - navdata['Toe'][idnan][bestephind].values # [datetime.datetime]
    # constants
    GM = 3986005.0E8 # universal gravational constant
    OeDOT = 7.2921151467E-5
    # Elements
    ecc = navdata['Eccentricity'][idnan][bestephind].values # Eccentricity
    Mk = navdata['M0'][idnan][bestephind].values + \
         t *(np.sqrt(GM / navdata['sqrtA'][idnan][bestephind].values**6) + 
             navdata['DeltaN'][idnan][bestephind].values)
    Ek = pyGnss.solveIter(Mk, ecc)
    Vk = np.arctan2(np.sqrt(1.0 - ecc**2) * np.sin(Ek), np.cos(Ek) - ecc)
    PhiK = Vk + navdata['omega'][idnan][bestephind].values
    # Perturbations
    delta_uk = navdata['Cuc'][idnan][bestephind].values * np.cos(2.0*PhiK) + \
              navdata['Cus'][idnan][bestephind].values * np.sin(2.0*PhiK)
    Uk = PhiK + delta_uk
    
    delta_rk = navdata['Crc'][idnan][bestephind].values * np.cos(2.0*PhiK) + \
               navdata['Crs'][idnan][bestephind].values * np.sin(2.0*PhiK)
    Rk = navdata['sqrtA'][idnan][bestephind].values**2 * (1.0 - ecc * np.cos(Ek)) + delta_rk
    
    delta_ik = navdata['Cic'][idnan][bestephind].values * np.cos(2.0*PhiK) + \
               navdata['Cis'][idnan][bestephind].values * np.sin(2.0*PhiK)
    Ik = navdata['Io'][idnan][bestephind].values + \
         navdata['IDOT'][idnan][bestephind].values * t + delta_ik
    
    #Compute the right ascension
    Omegak = navdata['Omega0'][idnan][bestephind].values + \
             (navdata['OmegaDot'][idnan][bestephind].values - OeDOT) * t - \
             (OeDOT * navdata['Toe'][idnan][bestephind].values)
             
    #X,Y coordinate corrections
    Xkprime = Rk * np.cos(Uk)
    Ykprime = Rk * np.sin(Uk)
    # ECEF XYZ
    X = Xkprime * np.cos(Omegak) - (Ykprime * np.sin(Omegak) * np.cos(Ik))
    Y = Xkprime * np.sin(Omegak) + (Ykprime * np.cos(Omegak) * np.cos(Ik))
    Z = Ykprime * np.sin(Ik)
    
    return np.array([X,Y,Z])

rtklibfn = os.getcwd() + '\\mah6.txt'

data = np.genfromtxt(rtklibfn, skip_header=1, delimiter='    ', dtype=str)
dt = np.array([datetime.strptime(t, '%Y/%m/%d %H:%M:%S.0') for t in data[:,0]])
sv = data[:,1]
az = data[:,2].astype(float)
el = data[:,3].astype(float)

prn = 'G23'
idsv = sv == prn

# PyGnss FOR MAH6
fnav = 'E:\\mahali\\2015-10-07\\brdc2800.15n'
fobs = 'E:\\mahali\\2015-10-07\\mah62800.15o.nc'
times = dt[idsv]

# Nav file
navdata = gr.load(fnav)
navdatatime = navdata.time.values
svclkbias = navdata.sel(sv=prn).SVclockBias.values
svBiasRange = svclkbias * 3e8
svrange = svBiasRange[0]

D = xarray.open_dataset(fobs, group='OBS', autoclose=True)
rx_xyz = D.position
C1 = D.sel(sv=prn)['C1'].values
obstimes64 = D.time.values
T = np.array([Timestamp(t).to_pydatetime() for t in obstimes64])
                        

aer = np.array(pyGnss.getSatellitePosition(rx_xyz, prn, times + timedelta(seconds=17), fnav, cs='aer',dtype='georinex'))
satxyz = np.array(pyGnss.getSatellitePosition(rx_xyz, prn, times + timedelta(seconds=17), fnav, cs='xyz',dtype='georinex'))
rho = np.sqrt( np.square(satxyz[0]-rx_xyz[0]) + np.square(satxyz[1]-rx_xyz[1]) + np.square(satxyz[2]-rx_xyz[2]) ) #- svrange
#New
satxyz = gpsSatPosition(fnav, times, sv=prn)
rec_lat, rec_lon, rec_alt = ecef2geodetic(rx_xyz[0], rx_xyz[1], rx_xyz[2])
A,E,R = ecef2aer(satxyz[0], satxyz[1], satxyz[2], rec_lat, rec_lon, rec_alt)

plt.figure(figsize=(10,7))
plt.title(prn)
plt.plot(dt[idsv], rho, 'b')
plt.plot(T, C1, 'r')

#plt.figure(figsize=(10,7))
#plt.plot(dt[idsv], el[idsv], '.b')
#plt.plot(dt[idsv], aer[1], 'r', lw=3)
#
plt.figure(figsize=(10,7))
plt.plot(dt[idsv], el[idsv]-aer[1], 'r')
plt.plot(dt[idsv], el[idsv]-E, 'k')
#
#plt.figure(figsize=(10,7))
#plt.plot(dt[idsv], el[idsv]-aer[1], 'r')

