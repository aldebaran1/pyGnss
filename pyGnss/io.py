# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 10:50:06 2018

@author: smrak
"""

import os
from glob import glob
import georinex as gr
import numpy as np
#from .gnssUtils import solveIter
from .pyGnss import getGpsTime, solveIter
from pymap3d import ecef2geodetic, ecef2aer

def getSatBias(fn: str = None):
    """
    Read satellite biases in TECu from JPL .XXi rinex file
    ftp://cddis.nasa.gov/gnss/products/ionex/
    Returns a dictionaty svbias, {'G05': -8.2} 
    """
    if os.path.isdir(fn):
        fn = glob(fn + '*.**i')
        if len(fn) >= 1:
            fn = fn[0] 
        else: 
            raise('No satbias files found in the folder')
            return 0
    i = 0
    svbias = {}
    with open(fn, 'r') as f:
        while 1:
            line = f.readline()
            if 'DIFFERENTIAL CODE BIASES' in line:
                while i < 32:
                    l = f.readline().split()[:2]
                    svbias['G'+l[0]] = float(l[1])
                    i += 1
                f.close()
                break
    return svbias

def getLeapSeconds(navfn):
    with open(navfn,'r') as f:
    #find end of header, which has non-constant length
        while True:
            line = f.readline()
            if 'LEAP SECONDS' in line:
                try:
                    leap_second = int(line.lstrip().rstrip().split()[0])
                    # Return
                    return leap_second
                except Exception as e:
                    print (e)
                    break
                break
            if 'END' in line:
                # Return null
                return 0
                break

def gpsSatPosition(fnav: str = None,
                   dt = None, 
                   sv: str = None, 
                   rx_position = None, 
                   coords: str = 'xyz'):
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
    gpstime = np.array([getGpsTime(t) for t in dt])
    t = gpstime - navdata['Toe'][idnan][bestephind].values # [datetime.datetime]
    # constants
    GM = 3986005.0E8 # universal gravational constant
    OeDOT = 7.2921151467E-5
    # Elements
    ecc = navdata['Eccentricity'][idnan][bestephind].values # Eccentricity
    Mk = navdata['M0'][idnan][bestephind].values + \
         t *(np.sqrt(GM / navdata['sqrtA'][idnan][bestephind].values**6) + 
             navdata['DeltaN'][idnan][bestephind].values)
    Ek = solveIter(Mk,ecc)
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
    
    if coords == 'xyz':
        return np.array([X,Y,Z])
    elif coords == 'aer':
        assert rx_position is not None
        rec_lat, rec_lon, rec_alt = ecef2geodetic(rx_position[0], rx_position[1], rx_position[2])
        A,E,R = ecef2aer(X, Y, Z, rec_lat, rec_lon, rec_alt)
        return np.array([A,E,R])
    
def singleRx(obs, nav, sv='G23', args=['L1','S1'], tlim=None,rawplot=False,
             sp=True,s4=False,polyfit=False,indicator=False,
             alt=300,el_mask=30,skip=20,porder=8,tec_ch=2,
             forder=5,fc=0.1,fs=1):
    
    def _plot(x,y,title=''):
        plt.figure()
        plt.title(title)
        plt.plot(x,y,'b')
    
    D = xarray.open_dataset(obs, group='OBS')
    rx_xyz = D.position
    leap_seconds = uf.getLeapSeconds(nav)
    obstimes64 = D.time.values
    times = np.array([Timestamp(t).to_pydatetime() for t in obstimes64]) - \
                        datetime.timedelta(seconds = leap_seconds)
    # define tlim ie. skip
    if tlim is not None:
        s = ((times>=tlim[0]) & (times<=tlim[1]))
    else:
        s = np.full(times.shape[0], True, dtype=bool)
        s[:skip] = False
    times = times[s]
    # Get satellite position
    aer = getSatellitePosition(rx_xyz, sv, times, nav, cs='aer',dtype='georinex')
    # Elevation mask
    idel = aer[1] >= el_mask
    times = times[idel]
    # times to output dict: Y
    Y = {'times': times}
    Y['rx_xyz'] = rx_xyz
    Y['az'] = aer[0][idel]
    Y['el'] = aer[1][idel]
    for arg in args:
        if not arg[1:] == 'TEC':
            X = D.sel(sv=sv)[arg].values[s][idel]
            # To dict
            Y[arg] = X
            # Indicators?
            if indicator:
                try:
                    if arg[0] == 'L' and arg[-1] != '5':
                        argi = arg + 'lli'
                    elif arg[0] == 'C' or arg[0] == 'P' or arg == 'L5':
                        argi = arg + 'ssi'
                    else:
                        argi = ''
                    Xi = D.sel(sv=sv)[argi].values[s:][idel]
                    Y[argi] = Xi
                except Exception as e:
                    print (e)
            if rawplot:
                _plot(times,X,arg)
        if arg[0] == 'L' or arg[0] == 'C':
            if arg == 'C1':
                X = X * f1 / c0 # To cycles
            Xd = phaseDetrend(X, order=porder)
                
            if polyfit:
                # To dict
                Y[arg+'polyfit'] = Xd
            if rawplot:
                if rawplot:
                    _plot(times,Xd,arg+'_detrend_poly')
            if sp:
                Xy = uf.hpf(Xd, order=forder,fc=fc,plot=False,fs=fs)
                # To dict
                Y['sp'] = Xy
                if rawplot:
                    _plot(times[100:],Xy[100:],arg+'_scint')
        if arg[0] == 'S':
            if s4:
                Xy = AmplitudeScintillationIndex(X,60)
                # To dict
                Y['s4'] = Xy
                if rawplot:
                    _plot(times,Xy,arg+'_S4')
                    
        if arg[1:] == 'TEC':
            if tec_ch == 2:
                C1 = D.sel(sv=sv)['C1'].values[s][idel]
                L1 = D.sel(sv=sv)['L1'].values[s][idel]
                C2 = D.sel(sv=sv)['P2'].values[s][idel]
                L2 = D.sel(sv=sv)['L2'].values[s][idel]
            elif tec_ch == 5:
                C1 = D.sel(sv=sv)['C1'].values[s][idel]
                L1 = D.sel(sv=sv)['L1'].values[s][idel]
                C2 = D.sel(sv=sv)['C5'].values[s][idel]
                L2 = D.sel(sv=sv)['L5'].values[s][idel]
            sTEC = getPhaseCorrTEC(L1, L2, C1, C2,channel=tec_ch)
            sTEC = sTEC - np.nanmin(sTEC)
            if arg == 'sTEC':
                # To dict
                Y[arg] = sTEC
                if rawplot:
                    _plot(times, sTEC,title=arg)
            elif arg == 'vTEC':
                vTEC = getVerticalTEC(sTEC,aer[1][idel],alt)
                # To dict
                Y[arg] = vTEC
                if rawplot:
                    _plot(times, vTEC, title=arg)
    # return dict with data/args
    return Y

def dataFromNC(fnc,fnav,sv,
               el_mask=30,tlim=None, 
               satpos=False,
               ipp=False,
               ipp_alt=None):
    leap_seconds = uf.getLeapSeconds(fnav)
    D = gr.load(fnc, useindicators=True).sel(sv=sv)
    if tlim is not None:
        if len(tlim) == 2:
            D = D.where(np.logical_and(D.time >= np.datetime64(tlim[0]), D.time <= np.datetime64(tlim[1])), drop=True)
    obstimes64 = D.time.values
    dt = np.array([Timestamp(t).to_pydatetime() for t in obstimes64]) - \
               datetime.timedelta(seconds = leap_seconds)
    rx_xyz = D.position
    aer = gpsSatPosition(fnav,dt,sv=sv, rx_position=rx_xyz, coords='aer')
    idel = (aer[1] >= el_mask)
    aer = aer[:, idel]
    dt = dt[idel]
    if satpos:
        D['az'] = aer[0]
        D['el'] = aer[1]
    if ipp:
        if ipp_alt is None:
            print ('Auto assigned altitude of the IPP: 250 km')
            ipp_alt = 250e3
        else:
            ipp_alt *= 1e3
        rec_lat, rec_lon, rec_alt = ecef2geodetic(rx_xyz[0], rx_xyz[1], rx_xyz[2])
        fm = np.sin(np.radians(aer[1]))
        r_new = ipp_alt / fm
        lla_vector = np.array(aer2geodetic(aer[0], aer[1], r_new, rec_lat, rec_lon, rec_alt))
        D['ipp_lon'] = lla_vector[1]
        D['ipp_lat'] = lla_vector[0]
        D['ipp_alt'] = ipp_alt
    return dt, D, idel