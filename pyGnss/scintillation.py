# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 14:18:35 2018

@author: smrak
"""

import numpy as np
from pyGnss import gnssUtils as gu
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

def getROTI(tec, length):
    """
    Sebastijan Mrak
    getROTI returns the rate of TEC Index calculated as the standard deviation 
    of the provided TEC on the moving window of the length 'length'. It returns 
    the ROTI as a numpy array data type.
    """
    rot = np.hstack((np.nan, np.diff(tec)))
    roti = sigmaTEC(rot, length)
    return roti
    

def phaseScintillationIndex(data, N):
    """
    Sebastijan Mrak
    GNSS Phase scintillation index for the interval of the length 'N' samples
    """
    y = np.nan * np.zeros(data.shape[0]-N)
    for i in range(data.shape[0] - N):
        if np.sum(np.isfinite(data[i:i+N])) > 5:
            y[i] = np.nanstd(data[i:i+N])
    return y
    
def AmplitudeScintillationIndex(data, N):
    """
    Sebastijan Mrak
    GNSS Amplitude scintillation index for the interval of the length 'N' samples
    """
    y = np.nan * np.zeros(data.shape[0])
    for i in range(data.shape[0] - N):
        if np.sum(np.isfinite(data[i:i+N])) > 2:
#            y[i] = np.sqrt( (np.nanmean(np.square(data[i:i+N]))-np.nanmean(data[i:i+N])**2)  / np.nanmean(data[i:i+N])**2 )
            y[i] = np.sqrt(np.nanvar(data[i:i+N])) / np.nanmean(data[i:i+N])
    return y

def sigmaTEC(x, N):
    idx = np.isnan(x)
    n2 = int(N/2)
    iterate = np.arange(n2, x.size-n2)
    y = np.nan * np.copy(x)
    for i in iterate:
        chunk = x[i-n2:i+n2]
        if np.sum(np.isfinite(chunk)) > N/4:
            y[i] = np.nanstd(chunk)
    y[idx] = np.nan
    return y

def s4(x, N):
    idx = np.isnan(x)
    n2 = int(N/2)
    iterate = np.arange(n2, x.size-n2)
    y = np.nan * np.copy(x)
    for i in iterate:
        chunk = x[i-n2:i+n2]
        if np.sum(np.isfinite(chunk)) > N/4:
            y[i] = np.nanstd(chunk) / np.nanmean(chunk)
    y[idx] = np.nan
    return y

def removeRipple(y, E = 5, L = 300, eps=False):
    std = np.nanstd(y[L:])
    envelope = gu.runningMax(y, N=60)
    std = np.nanmedian(envelope)
    e = E * std
    if np.where(abs(np.nan_to_num(y[:L])) >= e)[0].size > 0:
        if np.where(abs(np.nan_to_num(y[:L])) >= e)[0].size == 1:
            ex = np.where(abs(np.nan_to_num(y[:L])) >= e)[0].item() + 1
        else:
            ex = np.where(abs(np.nan_to_num(y[:L])) >= e)[0][-1] + 1
    else: 
        ex = -999
    
    if eps:
        return ex, e
    else:
        return ex

def mkrngs(y0, idf, gap_length=10, lim=0.05, min_length=None, max_length=None, 
            zero_mean=False, extend=0):
    gap = np.diff(np.where(idf)[0])
    i00 = np.where(idf)[0][0]
    i99 = np.where(idf)[0][-1]
    ixg = np.squeeze(np.argwhere(gap >= gap_length))
    LL = np.sort(np.hstack((ixg, ixg+1)))
    inner_limits = np.where(idf)[0][LL]
    limits = np.sort(np.hstack((i00,inner_limits,i99)))
    assert limits.size % 2 == 0
    ranges = limits.reshape(int(limits.size/2), 2)
    # Check for ranges vlidity: approx. zero mean
    if zero_mean:
        mask = []
        for i, r in enumerate(ranges):
            m_hat = np.nanmean(y0[r[0]:r[1]])
            if abs(m_hat) < lim: mask.append(i)
        if len(mask) > 0:
            mask = np.array(mask)
            ranges = ranges[mask]
    if min_length is not None:
        mask = np.squeeze(np.diff(ranges) > min_length)
        ranges = ranges[mask]
    if max_length is not None:
        mask = np.squeeze(np.diff(ranges) < max_length)
        ranges = ranges[mask]
    if len(ranges.shape) == 3 and ranges.shape[0] != 0: ranges = ranges[0]
    try:
        if extend > 0:
            start = ranges[:,0]
            ixstart = start > extend + 1
            ranges[ixstart,0] -= extend
            stop = ranges[:,1]
            ixstop = stop < (y0.size - extend - 1)
            ranges[ixstop, 1] += extend
    except:
        pass
    return ranges

def ranges(x, idf, min_gap=10, min_length=10*60, zero_mean=False):
    gap = np.diff(np.where(idf)[0])
    intervals = []
    if np.argwhere(gap >= min_gap).size > 0:
        intervals = mkrngs(x, idf, gap_length=min_gap, 
                            min_length=min_length, 
                            zero_mean=zero_mean)
    else:
        intervals = np.array([ [np.where(idf)[0][0], 
                             np.where(idf)[0][-1]+1] ])

    if len(intervals.shape) == 3: 
        try:
            intervals = intervals[0]
        except: 
            intervals = np.array([])
    
    return intervals

def partialProcess(dt,r, x, fs=1, fc=0.1, hpf_order=6,
                    remove_outliers = False,
                    plot_ripple = False,
                    plot_outlier = False):
    
    idf = np.isfinite(x)
    # If there are NaNs in the interval, do a cubic spline. 
    # Max gap is 10 seconds set by the "make ranges routine"
    # 1. dTEC Split
    if np.sum(np.isnan(x)) > 0:
        x0 = np.where(idf)[0]
        x1 = np.arange(x.size)
        CSp = CubicSpline(x0, x[idf])
        x_cont = CSp(x1)
    else:
        x_cont = np.copy(x)
    # 2. Tec/snr scintillation (high-pass) filtering!
    tec_hpf = gu.hpf(x_cont, fs=fs, order=hpf_order, fc=fc)
    tec_hpf[~idf] = np.nan
    
    # 3. Remove initial ripple on the scintillation time-series
    sT_exit, eps = removeRipple(tec_hpf, E=1.5, L=300, eps = True)
    if plot_ripple:
        plt.figure()
        plt.plot(dt[r[0]:r[1]], tec_hpf, 'b')
        plt.plot([dt[r[0]], dt[r[1]]], [eps, eps], '--r')
        if sT_exit != -999: plt.plot(dt[r[0]:r[1]][:sT_exit], tec_hpf[:sT_exit], 'xr')
    if sT_exit != -999: 
        tec_hpf[:sT_exit] = np.nan
    # 4. Outlier detection and removal. Still on the scintillation time-series.
    # 4.1 TEC Scintillation
    if remove_outliers:
        envelope = gu.runningMax(abs(tec_hpf), N = 10)
        median_envelope = gu.runningMedian(envelope, N = 120)
        outlier_margin = median_envelope + 5 * np.nanstd(tec_hpf)
        idoutlier = np.nan_to_num(abs(tec_hpf)) > outlier_margin
        tec_hpf[idoutlier] = np.nan
        if np.nansum(idoutlier) > 0:
            if plot_outlier:
                plt.figure(figsize = [8,5])
                plt.plot(dt[r[0] : r[1]], tec_hpf, 'b', label = '$\delta TEC_{0.1 Hz}$')
                plt.plot(dt[r[0] : r[1]], median_envelope, 'g', label = 'env = <$\widehat{\delta TEC}>|_{10s}$')
                plt.plot(dt[r[0] : r[1]], outlier_margin, '--r', label = '$\epsilon$ = env + 4$\cdot \sigma(\delta TEC)|_{60s}$')
                plt.plot(dt[r[0] : r[1]], -outlier_margin, '--r')
                plt.plot(dt[r[0] : r[1]][~idoutlier], tec_hpf[~idoutlier], 'xr')
                plt.ylabel('$\delta$ TEC [TECu]')
                plt.xlabel('Time [UTC]')
                plt.grid(axis='both')
                plt.legend()
    return tec_hpf
