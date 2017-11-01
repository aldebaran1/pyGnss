#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 14:42:15 2017

@author: Sebastijan Mrak <smrak@gmail.com>
"""

from scipy import signal
import numpy as np
import os
import glob
import datetime
import pandas
from pandas import read_hdf
import yaml
import h5py
from mpl_toolkits.basemap import Basemap
from gsit import pyGps
from pyGnss import pyGnss 
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def butter_hpf(highcut, fs, order):
    """
    Sebastijan Mrak
    Design the Butterwoth response highpass filter with N-th order and 
    3db cutoff frequency 'highcut' in Hz.
    Output are the poles 'b' and zeroes 'a' of the filter
    """
    nyq = 0.5 * fs
    high = highcut / nyq
    b, a = signal.butter(order, high, btype='highpass')
    w, h = signal.freqz(b, a, worN=1000)
    return b, a

def butter_lpf(fc, fs, order):
    """
    Sebastijan Mrak
    Design the Butterwoth response highpass filter with N-th order and 
    3db cutoff frequency 'highcut' in Hz.
    Output are the poles 'b' and zeroes 'a' of the filter
    """
    nyq = 0.5 * fs
    high = fc / nyq
    b, a = signal.butter(order, high, btype='lowpass', analog=False)
    w, h = signal.freqz(b, a, worN=1000)
    return b, a


def bpf(y, lowcut, highcut, fs=1, order=5, plot=False):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    mid = ((highcut - lowcut)/3)*2 / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    w, h = signal.freqz(b, a, worN=1000)
    gd = -np.diff(np.unwrap(np.angle(h)))/np.diff(w)
    idx = abs(w-high).argmin()
    y_filt = signal.lfilter(b, a, y)
    if plot:
        plt.figure()
        plt.semilogx(w, 20*np.log10(np.abs(h)))
        plt.plot(w[idx], 20*np.log10(np.abs(h[idx])), 'xr')
        plt.ylim([-60,5])
        plt.title('Magnitude-normalized Butterworth filter frequency response')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Amplitude [dB]')
        plt.grid(which='both', axis='both')
        ################################################
    
        plt.figure()
        plt.semilogx(w[1:], gd)
        plt.plot(w[idx], gd[idx], 'xr')
        plt.title('LPF group delay')
        plt.xlabel('Frequency [radians / second]')
        plt.ylabel('Group delay [samples]')
        plt.margins(0, 0.1)
        plt.grid(which='both', axis='both')
    return y_filt , gd[idx]

def hpf(y, fc=0.1, order=5, fs=1,plot=False, group_delay=False):
    """
    Sebastijan Mrak
    Filter the input data 'y' with desired HP filter.  
    """
    b, a = butter_hpf(fc, fs, order)
    y_filt = signal.lfilter(b, a, y)
    w, h = signal.freqz(b, a, worN=1000)
    gd = -np.diff(np.unwrap(np.angle(h)))/np.diff(w)
    print ('Group delay of the filter is '+ str(gd[-1])+' samples.')
    if plot:
        plt.figure()
        plt.semilogx(w, 20*np.log10(np.abs(h)))
        plt.ylim([-60,5])
        plt.title('Magnitude-normalized Butterworth filter frequency response')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Amplitude [dB]')
        plt.grid(which='both', axis='both')
        ################################################
    
        plt.figure()
        plt.semilogx(w[1:], gd)
        plt.title('LPF group delay')
        plt.xlabel('Frequency [radians / second]')
        plt.ylabel('Group delay [samples]')
        plt.margins(0, 0.1)
        plt.grid(which='both', axis='both')
    if group_delay:
        return y_filt, gd[-1]
    else:
        return y_filt

def lpf(y, fc=0.1, order=5, fs=1, plot=False, group_delay=False):
    b, a = butter_lpf(fc, fs, order)
    y_filt = signal.lfilter(b, a, y)
    w, h = signal.freqz(b, a, worN=1000)
    gd = -np.diff(np.unwrap(np.angle(h)))/np.diff(w)
    print ('Group delay of the filter is '+ str(gd[1])+' samples.')
    if plot:
        plt.figure()
        plt.semilogx(w, np.log10(np.abs(h)))
        plt.ylim([-60,5])

        plt.title('Magnitude-normalized butterworth filter frequency response')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Amplitude [dB]')
        plt.grid(which='both', axis='both')
        ################################################

        plt.figure()
        plt.semilogx(w[1:], gd)
        plt.title('LPF group delay')
        plt.xlabel('Frequency [radians / second]')
        plt.ylabel('Group delay [samples]')
        plt.margins(0, 0.1)
        plt.grid(which='both', axis='both')
    if group_delay:
        return y_filt, gd[1]
    else:
        return y_filt
#------------------------------------------------------------------------------#
def interpolateLatLon(lat, lon, order=3, resolution=0.1):
    z = np.polyfit(lon, lat, order)
    f = np.poly1d(z)
    x_new = np.arange(min(lon), max(lon), resolution)
    y_new = f(x_new)
    return y_new, x_new

def interpolateTEC(t,tec,order=20,maxjump=2):
    idf = np.where(np.isfinite(tec))[0]
    if np.diff(tec[idf]).max() <= maxjump:
        x = np.arange(0,tec.shape[0])
        z = np.polyfit(x[np.isfinite(tec)], tec[np.isfinite(tec)], order)
        f = np.poly1d(z)
        y = f(x)
        mask = np.where(np.isnan(tec))[0]
        tec[mask] = y[mask]
#    print (np.diff(mask))
    return tec

#------------------------------------------------------------------------------#
def polynom(y, order=3):
    x = range(y.shape[0])
    z = np.polyfit(x, y, order)
    f = np.poly1d(z)
    y_new = f(x)
    return y_new
#------------------------------------------------------------------------------#
def correctSampling(t, y, fs=1):
    ts = pyGps.datetime2posix(t)
    td = np.diff(ts)
    idt = np.where(td != fs)[0]
    if idt.shape[0] > 0:
        while True:
            td = np.diff(ts)
            idt = np.where(td != fs)[0]
            if idt.shape[0] == 0:
                break
            ts = np.insert(ts, idt[0]+1, ts[idt[0]]+fs)
            y = np.insert(y, idt[0]+1, np.NaN)
    return ts, y
#------------------------------------------------------------------------------#
def returnSlope(t,y, fs=5, interval=5):
    skip = int(60/fs) * interval
    t_new = t[::skip]
    slope = np.diff(y[::skip])
    return t_new[:-1], slope
#------------------------------------------------------------------------------#
def getIntervals(y, maxgap=1, maxjump=0.5):
    r = np.array(range(len(y)))
    idx = np.isfinite(y)
    r = r[idx]
    intervals=[]
    if len(r)==0:
        return intervals
    beginning = 0
    last = 0
    for i in range(len(r)):
        if (r[i]-r[last] > maxgap) or (abs(y[i] - y[last]) > maxjump):
            intervals.append((r[beginning],r[last]))
            beginning=i
        last=i
        if r[i]==r[-1]:
            intervals.append([r[beginning],r[last]])
            break
    return intervals
#------------------------------------------------------------------------------#
def getWegihtedPolyfit(ix, t, tec, Tdelta=10, Ts=30, interval_mode=2, overlap = 5):
    ix_min = int(20 * 60 / Ts) #minutes normailized to samples

    if ix > ix_min and ix < tec.shape[0] - ix_min:
        
#        Tdelta = 10 #minutes
        # find the interval
        t1 = t[ix] - datetime.timedelta(minutes=Tdelta)
        t2 = t[ix] + datetime.timedelta(minutes=Tdelta)
        t1_ix = abs(t-t1).argmin()
        t2_ix = abs(t-t2).argmin()
        
        if interval_mode == 3:
            intervals = [(0, t1_ix+5), (t1_ix-5, t2_ix+5), (t2_ix-5, len(t))]#, (t2_ix, len(t))]
        else:
            intervals = [(0, t1_ix), (t1_ix, len(t))]
            
        polyfit = np.nan*np.ones(tec.shape[0])
        c=0
        for lst in intervals:
            # Given a length of the interval decide what order of polynom it should be
            polynom_order = getPolynomOrder(lst, Ts)
            # Init am temp TEC subarray array and get continous intervals 
            if c == 0:
                tmp = tec[lst[0]:lst[1]+overlap]
                c+=1
            else:
                tmp = tec[lst[0]-overlap:lst[1]]
                c+=1
            LST = getIntervals(tmp,maxgap=1, maxjump=1)

            # Init an empty polynom array of a lenght of the TEC array
            pp = np.nan*np.ones(tmp.shape[0])
            for l in LST:
                # Ignore intervals with a length of less than 10 successive samples
                if l[1]-l[0] > 10:
                    polynom_order = getPolynomOrder(l, Ts)
                    # Get a polyfit for a given interval
                    pp[l[0]:l[1]] = polynom(tmp[l[0]:l[1]], order=polynom_order)
            # Save polynoms per interval into common polynom array
            if c == 1:
                polyfit[lst[0]:lst[1]+overlap] = pp
            else:
                polyfit[lst[0]-overlap:lst[0]] = np.add(polyfit[lst[0]-overlap:lst[0]], pp[0:overlap])/2
                polyfit[lst[0]:lst[1]] = pp[overlap:]
            
    else:
        intervals = getIntervals(tec, maxgap=1, maxjump=1)
        pp = np.nan*np.ones(tec.shape[0])
        for lst in intervals:
            if lst[1]-lst[0] > 10:
                polynom_order = getPolynomOrder(lst, Ts)
                pp[lst[0]:lst[1]] = polynom(tec[lst[0]:lst[1]], order=polynom_order)
        polyfit = pp
    polyfit[:10] = np.nan
    polyfit[-10:] = np.nan
    return polyfit

#------------------------------------------------------------------------------#
def _alignTimes(tlist, teclist, polylist, residuallist, fs):
    tmin = []
    tmax = []
    for i in range(len(tlist)):
        tmin.append(tlist[i].min())
        tmax.append(tlist[i].max())
    tstart = max(pyGps.datetime2posix(tmin))
    tend = min(pyGps.datetime2posix(tmax))
    
    t = []
    tec2 = []
    poly2 = []
    res2 = []
    for i in range(len(teclist)):
        tt, tec1 = correctSampling(tlist[i], teclist[i], fs=fs)
        tt, poly1 = correctSampling(tlist[i], polylist[i], fs=fs)
        tt, res1 = correctSampling(tlist[i], residuallist[i], fs=fs)
        tt = np.array(tt)
        idt = np.where((tt>=tstart) & (tt<=tend))[0]
        t.append(tt[idt])
        tec2.append(tec1[idt])
        poly2.append(poly1[idt])
        res2.append(res1[idt])
    return t, tec2, poly2, res2
#------------------------------------------------------------------------------#
def getRxList(folder, sufix):
    filestr = os.path.join(folder,sufix)
    flist = sorted(glob.glob(filestr))
    rx = []
    for f in flist:
        head, tail = os.path.split(f)
        rx.append(tail[0:4])
    return rx
#------------------------------------------------------------------------------#
def createTimeArray(timelim):
    ts = pyGps.datetime2posix(timelim)
    t = range(int(ts[0]), int(ts[1])+1)
    return np.array(t)
#------------------------------------------------------------------------------#
def returnTotalityPath(width=False):
    data = h5py.File('totality.h5', 'r')
    time = np.array(data['path/time'])
    center_lat = np.array(data['path/center_lat'])
    center_lon = np.array(data['path/center_lon'])
    north_lat = np.array(data['path/north_lat'])
    north_lon = np.array(data['path/north_lon'])
    south_lat = np.array(data['path/south_lat'])
    south_lon = np.array(data['path/south_lon'])
    
    if not width:
        return time, center_lat, center_lon
    else:
        return time, [center_lat, north_lat, south_lat], [center_lon, north_lon, south_lon]
#------------------------------------------------------------------------------#
def EuclidDistance(lon1, lat1, lon2, lat2):
    yy = (lat1-lat2) * 111 # in km for LAT diff
    xx = (lon1-lon2) * 111 / np.cos(np.radians(lat2)) #In km with lat angle correction
    distance = np.sqrt( pow(xx,2) + pow(yy,2))
    if lat2 > lat1:
        return distance, xx, yy
    else:
        return distance, xx, yy
#------------------------------------------------------------------------------#
def getToatlityTouch(totality, los):
    euclid_dist = []
    errX = []
    Ix = []
    Tt = totality[0]
    Tipp = los[0]
    for i in range(len(Tt)):
        Ts = Tt[i]
        j = abs(Ts - Tipp).argmin()
        ed, ex, ey = EuclidDistance(totality[2][i], totality[1][i], los[2][j], los[1][j])
        euclid_dist.append(ed)
        errX.append(ex)
        Ix.append(j)
    ixe = abs(np.array(euclid_dist)).argmin()
    ix = Ix[ixe-1]
    return ix, np.array(euclid_dist), errX
#------------------------------------------------------------------------------#
def getPolynomOrder(interval, Ts, weights=[2,3,7,10]):
    lst_len = (interval[1] - interval[0]) * Ts / 60 # In minutes
    
    if lst_len < 20:
        polynom_order = weights[0]
    elif lst_len >= 20 and lst_len < 50:
        polynom_order = weights[1]
    elif lst_len >= 50 and lst_len < 100:
        polynom_order = weights[2]
    elif lst_len >= 100:
        polynom_order = weights[3]
    return polynom_order
#------------------------------------------------------------------------------#
def resBPFFilter(t, y, Ts=30):
    bpf_res = np.nan*np.ones(y.shape[0])
    idx = np.where(np.isfinite(y))[0]
    tmp_res, gd = bpf(y[idx], 0.0001, 0.002, fs=1/30, order=5)
    bpf_res[idx] = tmp_res
    bpf_td = np.array(t) - datetime.timedelta(seconds=9*Ts)
    
    return bpf_td, bpf_res
#------------------------------------------------------------------------------#

def getRxListCoordinates(filename='corsrxlist.h5', FDIR = '/home/smrak/Documents/eclipse/'):

    with (h5py.File(FDIR+filename, 'r')) as f:
        names = f['/data/rx']
        data = f['/data/table']
        rx = [rx[0].decode('utf-8') for rx in names]
        lat = [float(l[0]) for l in data]
        lon = [float(l[1]) for l in data]
    
    rxpos = np.array((lat, lon))
    return rx, rxpos

def rxFilter(rxlist, rxl, latlim=[], lonlim=[]):
    if isinstance(rxlist, list):
        rxlist = np.array(rxlist)
    elif isinstance(rxlist, np.ndarray):
        pass
    else:
        print ('Enter Rxlist as numpy array or a list')
        return
    
    if (len(latlim) == 0) and (len(lonlim) > 0):
        idx = np.where( (rxl[1] >= lonlim[0]) & (rxl[1] <= lonlim[1]) )[0]
        rxl_f = rxl[:,idx]
        rx = rxlist[idx]
    elif (len(latlim) > 0) and (len(lonlim) == 0):
        idx = np.where( (rxl[0] >= latlim[0]) & (rxl[0] <= latlim[1]) )[0]
        rxl_f = rxl[:,idx]
        rx = rxlist[idx]
    elif (len(latlim) > 0) and (len(lonlim) > 0):
        idx = np.where( (rxl[0] >= latlim[0]) & (rxl[0] <= latlim[1]) &
                        (rxl[1] >= lonlim[0]) & (rxl[1] <= lonlim[1]))[0]
        rxl_f = rxl[:,idx]
        rx = rxlist[idx]
    else:
        print ('Enter filtering limits!')
        return [], []
    
    return rx, rxl_f
################################################################################
def getPhaseCorrTECGLONASS(L1,L2,P1,P2,fN):
    f1 = (1602 + fN*0.5625) * 1000000
    f2 = (1246 + fN*0.4375) * 1000000
    c0 = 3E8
    range_tec = ((f1**2 * f2**2) / (f1**2 - f2**2)) * (P2 - P1) / 40.3 / pow(10, 16)
    phase_tec = ((f1**2 * f2**2) / (f1**2 - f2**2)) * (c0/40.3) * (L1 / f1 - L2 / f2) / pow(10, 16)
    
    tec_difference = np.array(sorted(phase_tec-range_tec))
            
    tec_difference = tec_difference[np.isfinite(tec_difference)]
    median_difference = tec_difference[int(len(tec_difference)/2)]
    tec = phase_tec - median_difference
    
    return tec

def RxbiasEstimator(RxBtime=0, TECref=0, t=0, sTEC=0, sat_bias=0, el=[], IPPalt=350):
    Fm = pyGnss.getMappingFunction(el, IPPalt)
    
    idt = abs(t - RxBtime).argmin()
    Fmtt = Fm[idt]
    sTECt = sTEC[idt]
    
    rx_bias = TECref/Fmtt - sTECt - sat_bias
    
    return rx_bias    


# ---------------------------------------------------------------------------- #
def returnTEC(data, sv, navfile, yamlfile, timelim=None, el_mask=30, leap_seconds=18, fN=0,
              alt=300, sattype='G', el=False, lla=False, vertical=False, svbias=False, raw=False,
              rxbias=False, RxB=[]):
    obstimes = np.array((data.major_axis))
    obstimes = pandas.to_datetime(obstimes) - datetime.timedelta(seconds=leap_seconds)
    stream = yaml.load(open(yamlfile, 'r'))
    rx_xyz = stream.get('APPROX POSITION XYZ')
    if (timelim is not None):
        idt = np.where( (obstimes>timelim[0]) & (obstimes<timelim[1]) ) [0]
        t = obstimes[idt]
        L1 = np.array(data['L1', sv, :, 'data'])
        L2 = np.array(data['L2', sv, :, 'data'])
        C1 = np.array(data['C1', sv, :, 'data'])
        C2 = np.array(data['P2', sv, :, 'data'])
        L1 = L1[idt]
        L2 = L2[idt]
        C1 = C1[idt]
        C2 = C2[idt]
        if sattype == 'R':
            aer = pyGnss.getIonosphericPiercingPoints(rx_xyz, sv-32, t, alt, navfile, cs='aer', sattype=sattype)
            llt = pyGnss.getIonosphericPiercingPoints(rx_xyz, sv-32, t, alt, navfile, cs='wsg84', sattype=sattype)
        else:
            aer = pyGnss.getIonosphericPiercingPoints(rx_xyz, sv, t, alt, navfile, cs='aer', sattype='G')
            llt = pyGnss.getIonosphericPiercingPoints(rx_xyz, sv, t, alt, navfile, cs='wsg84', sattype='G')
##        filter by Elevation angle

        idel = np.where((aer[1] > el_mask))[0]     
        
        if vertical == False:
            tec = pyGps.getPhaseCorrTEC(L1[idel],L2[idel], C1[idel], C2[idel])
        else:
            tec = pyGps.getPhaseCorrTEC(L1[idel],L2[idel], C1[idel], C2[idel])
            if svbias and not rxbias:
                bstream = yaml.load(open('/media/smrak/Eclipse2017/Eclipse/jplg2330.yaml', 'r'))
                sat_bias = float(bstream.get(sv))
                tec = pyGnss.getVerticalTEC(tec+sat_bias, aer[1][idel], alt)
            if rxbias and svbias == False:
                tec = pyGps.getPhaseCorrTEC(L1[idel],L2[idel], C1[idel], C2[idel])
                bstream = yaml.load(open('/media/smrak/Eclipse2017/Eclipse/jplg2330.yaml', 'r'))
                sat_bias = float(bstream.get(sv))
                rx_bias = RxbiasEstimator(RxB[0], RxB[1], t[idel], tec, sat_bias, el=aer[1][idel], IPPalt=350)
                print (rx_bias)
                tec = pyGnss.getVerticalTEC(tec+sat_bias+rx_bias, aer[1][idel], alt)
            if svbias == False and rxbias == False:
                tec = pyGnss.getPhaseCorrTEC(L1[idel],L2[idel],C1[idel],C2[idel])
                tec = pyGnss.getVerticalTEC(tec, aer[1][idel], alt)
                
        t = t[idel]
    else:
        L1 = np.array(data['L1', sv, :, 'data'])
        L2 = np.array(data['L2', sv, :, 'data'])
        C1 = np.array(data['C1', sv, :, 'data'])
        C2 = np.array(data['P2', sv, :, 'data']) 
        t = obstimes
        if sattype == 'R':
            aer = pyGnss.getIonosphericPiercingPoints(rx_xyz, sv-32, t, alt, navfile, cs='aer', sattype=sattype)
            llt = pyGnss.getIonosphericPiercingPoints(rx_xyz, sv-32, t, alt, navfile, cs='wsg84', sattype=sattype)
        else:
            aer = pyGnss.getIonosphericPiercingPoints(rx_xyz, sv, t, alt, navfile, cs='aer', sattype='G')
            llt = pyGnss.getIonosphericPiercingPoints(rx_xyz, sv, t, alt, navfile, cs='wsg84', sattype='G')
##        filter by Elevation angle
        idel = np.where((aer[1] > el_mask))[0]

        if vertical == False:
            tec = pyGps.getPhaseCorrTEC(L1[idel],L2[idel], C1[idel], C2[idel])
        else:
            tec = pyGps.getPhaseCorrTEC(L1[idel],L2[idel], C1[idel], C2[idel])
            if svbias and not rxbias:
                bstream = yaml.load(open('/media/smrak/Eclipse2017/Eclipse/jplg2330.yaml', 'r'))
                sat_bias = float(bstream.get(sv))
                tec = pyGnss.getVerticalTEC(tec+sat_bias, aer[1][idel], alt)
            if rxbias and svbias == False:
                tec = pyGps.getPhaseCorrTEC(L1[idel],L2[idel], C1[idel], C2[idel])
                bstream = yaml.load(open('/media/smrak/Eclipse2017/Eclipse/jplg2330.yaml', 'r'))
                sat_bias = float(bstream.get(sv))
                rx_bias = RxbiasEstimator(RxB[0], RxB[1], t[idel], tec, sat_bias, el=aer[1][idel], IPPalt=350)
                print (rx_bias)
                tec = pyGnss.getVerticalTEC(tec+sat_bias+rx_bias, aer[1][idel], alt)
            if svbias == False and rxbias == False:
                tec = pyGnss.getPhaseCorrTEC(L1[idel],L2[idel],C1[idel],C2[idel])
                tec = pyGnss.getVerticalTEC(tec, aer[1][idel], alt)
        t = t[idel]
        
    #Return time, TEC and respective third argumnet masked by an elavation mask angle
    if el==False and lla==False and raw==False:
        return t, tec
    elif el==True and lla==False and raw==False:
        return t, tec, aer[:,idel]
    elif el==False and lla ==True and raw==False:
        return t, tec, llt[:,idel]
    elif el==False and lla ==True and raw==True:
        return t, tec, llt[:,idel], [L1[idel],L2[idel],C1[idel],C2[idel]]
    elif el==True and lla ==False and raw==True:
        return t, tec, aer[:,idel], [L1[idel],L2[idel],C1[idel],C2[idel]]
    else:
        return t, tec, aer[:,idel]
        
            
    
################################################################################
def getResiduals(rxlist=['mobu'],decimate='_30',sv=2,day=233,Ts=30,el_mask=20,
                 polynom_order=12, timelim=None, latlon=False):

    if timelim is None:
        timelim = [datetime.datetime.strptime('2017 '+str(day)+' 15 0 0', '%Y %j %H %M %S'), 
                   datetime.datetime.strptime('2017 '+str(day)+' 21 0 0', '%Y %j %H %M %S')]

    res = []
    time = []
    wsg = []
    for rx in rxlist:
        folder = '/media/smrak/Eclipse2017/Eclipse/cors/all/'+str(day)+'/'
        navfile = '/media/smrak/Eclipse2017/Eclipse/nav/jplm'+str(day)+'0.17n'
        hdffile = folder + rx + str(day)+'0'+decimate+'.h5'
        yamlfile = folder + rx + str(day)+'0'+decimate+'.yaml'
        
        
        try:
            data = read_hdf(hdffile)
            t, tec, lla, raws = returnTEC(data, sv=sv, navfile=navfile, yamlfile=yamlfile, 
                                             timelim=timelim, el_mask=el_mask, lla=True, 
                                             svbias=True, vertical=True, raw=True)
            
            ix, intervals = getIntervals(tec, maxgap=1, maxjump=1)
            p = np.nan*np.ones(tec.shape[0])
            for lst in intervals:
                p[lst[0]:lst[1]] = polynom(tec[lst[0]:lst[1]], order=polynom_order)
            p[0:10] = np.nan
            p[-10:] = np.nan
            z = tec - p
            res.append(z)
            time.append(t)
            wsg.append(lla)
        except:
            print (rx)
    if latlon:
        return time, res, wsg
    else:
        return time, res
################################################################################

################################################################################
def _plotLOS(tlist, teclist, polylist, residuallist, rx='', sv=0, save=False,
             fig_path=None,
             pltlim = [datetime.datetime(2017,8,21,16,0,0), datetime.datetime(2017,8,21,21,0,0)]):
    
    fig = plt.figure(figsize=(12,8))
    tdt = [datetime.datetime.utcfromtimestamp(i) for i in tlist[1]]
    formatter = mdates.DateFormatter('%H:%M')
    
    ax1 = fig.add_subplot(411)
    ax2 = fig.add_subplot(412, sharex=ax1)
    ax3 = fig.add_subplot(413, sharex=ax1)
    ax4 = fig.add_subplot(414, sharex=ax1)
    
    ax1.plot(tdt, teclist[0], 'b')
    ax1.plot(tdt, polylist[0], 'r')
    ax1.plot(tdt, teclist[1], 'b')
    ax1.plot(tdt, polylist[1], 'r')
    
    ax2.plot(tdt, residuallist[0], 'g')
    ax2.plot(tdt, residuallist[1], 'm')
    ax2.plot( [tdt[0], tdt[-1]], [0,0], '--k')
    
    ax3.plot(tdt, polylist[0], 'r')
    ax3.plot(tdt, teclist[1], 'b')

    ax4.plot(tdt, teclist[1]-polylist[0], 'm')
    
    plt.setp(ax1.get_xticklabels(), visible=False) 
    plt.setp(ax2.get_xticklabels(), visible=False) 
    plt.setp(ax3.get_xticklabels(), visible=False) 
    
    ax4.grid(axis='y')
    ax2.grid(axis='y')
    ax1.set_title(rx+' - sv: '+str(sv))
    ax1.set_ylabel('vTEC')
    ax2.set_ylabel('residuals')
    ax3.set_ylabel('vTEC')
    ax4.set_ylabel('diff TEC')
    ax4.set_xlabel('UTC')
    ax3.set_xlim(pltlim)
    ax4.xaxis.set(major_formatter=formatter)
    ax2.set_ylim([-0.5, 0.5])
    plt.tight_layout()
    if save == True:
        if fig_path is None:
            plt.savefig('/media/smrak/Eclipse2017/Eclipse/plots/cors/run4/'+rx+'_'+str(sv)+'.png', dpi=300)
        else:
            plt.savefig(fig_path+rx+'_'+str(sv)+'.png', dpi=300)
        plt.close(fig)
    else:
        plt.show()
#------------------------------------------------------------------------------#
def _plotEclipseMap(filepath='totality.h5'):
    data = h5py.File(filepath, 'r')
    center_lat = np.array(data['path/center_lat'])
    center_lon = np.array(data['path/center_lon'])
    north_lat = np.array(data['path/north_lat'])
    north_lon = np.array(data['path/north_lon'])
    south_lat = np.array(data['path/south_lat'])
    south_lon = np.array(data['path/south_lon'])
    
    (fig,ax) = plt.subplots(1,1,figsize=(16,12),facecolor='w')
    latlim2 = [33, 38]
    lonlim2 = [-95, -75]
    m = Basemap(projection='merc',
    llcrnrlat=latlim2[0],urcrnrlat=latlim2[1],\
    llcrnrlon=lonlim2[0],urcrnrlon=lonlim2[1],\
    resolution='c')
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    
    X,Y = m(center_lon, center_lat)
    X1,Y1 = m(north_lon, north_lat)
    X2,Y2 = m(south_lon, south_lat)
    m.plot(X,Y, c='r')
    m.plot(X1,Y1, c='b')
    m.plot(X2,Y2, c='b')
    plt.show()
#------------------------------------------------------------------------------#
def _plotDetrending(t, z, sv='',rx=None,order=[],save=False,fig_path=None, 
                    pltlim = [datetime.datetime(2017,8,21,17,0,0), datetime.datetime(2017,8,21,21,0,0)]):
    
    fig = plt.figure(figsize=(12,8))
    formatter = mdates.DateFormatter('%H:%M')
    axs = []
    for i in range(len(z)):
        ax = fig.add_subplot(len(z),1,i+1)
        axs.append(ax)
        plt.setp(ax.get_xticklabels(), visible=False) 
        plt.plot(t, z[i], label=order[i])
        plt.plot( [pltlim[0], pltlim[-1]], [0,0], '--k')
        plt.ylim([-0.8, 0.6])
        plt.legend()
        ax.set_ylabel('residuals')
        ax.grid(axis='y')
        ax.set_xlim(pltlim)
    ax.xaxis.set(major_formatter=formatter)
    plt.setp(ax.get_xticklabels(), visible=True)
    ax.set_xlabel('UTC')
    axs[0].set_title(rx+' - sv: '+str(sv))
    plt.tight_layout()
    
    if save == True:
        if fig_path is None:
            plt.savefig('/media/smrak/Eclipse2017/Eclipse/plots/cors/detrending/orders_'+rx+'_'+str(sv)+'.png', dpi=300)
        else:
            plt.savefig(fig_path+rx+'_'+str(sv)+'.png', dpi=300)
        plt.close(fig)
    else:
        plt.show()
        plt.close(fig)
#------------------------------------------------------------------------------#
def _plotResidualsTEC(t,vtec,stec,residuallist1,residuallist2,sv='',rx=None,order=[],save=False,fig_path=None,
                   pltlim = [datetime.datetime(2017,8,21,17,0,0), datetime.datetime(2017,8,21,21,0,0)]):
    fig = plt.figure(figsize=(12,8))
    formatter = mdates.DateFormatter('%H:%M')
    
    ax1 = fig.add_subplot(411)
    ax2 = fig.add_subplot(412, sharex=ax1)
    ax3 = fig.add_subplot(413, sharex=ax1)
    ax4 = fig.add_subplot(414, sharex=ax1)
    
    ax1.plot(t, stec, 'b')
    
    for i in range(len(residuallist1)):
        ax2.plot(t,residuallist1[i], label = order[i])
    ax2.plot( [t[0], t[-1]], [0,0], '--k')
    ax2.grid(axis='y')
    
    ax3.plot(t, vtec, 'b')
    
    for i in range(len(residuallist1)):
        ax4.plot(t,residuallist2[i], label = order[i])
    ax4.plot( [t[0], t[-1]], [0,0], '--k')
    ax4.grid(axis='y')
    
    ax1.set_title(rx+' - sv: '+str(sv))
    ax1.set_ylabel('sTEC')
    ax2.set_ylabel('residuals')
    ax3.set_ylabel('vTEC')
    ax4.set_ylabel('residuals')
    ax4.set_xlabel('UTC')
    plt.legend()
    ax4.xaxis.set(major_formatter=formatter)
    plt.setp(ax1.get_xticklabels(), visible=False) 
    plt.setp(ax2.get_xticklabels(), visible=False) 
    plt.setp(ax3.get_xticklabels(), visible=False) 
    ax4.set_xlim(pltlim)
    
    plt.tight_layout()
    if save == True:
        if fig_path is None:
            plt.savefig('/media/smrak/Eclipse2017/Eclipse/plots/cors/detrending/'+rx+'_'+str(sv)+'.png', dpi=300)
        else:
            plt.savefig(fig_path+rx+'_'+str(sv)+'.png', dpi=300)
        plt.close(fig)
    else:
        plt.show()
        plt.close(fig)
#------------------------------------------------------------------------------#
def wavesPlot(time, res, dist_coeff=0.5):
    plt.figure()
    for i in range(len(res)):
        plt.plot(time[i], res[i]+i*0.15, 'b')
