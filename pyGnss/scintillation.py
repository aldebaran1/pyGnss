# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 14:18:35 2018

@author: smrak
"""

import numpy as np

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
    
def AmplitudeScintillationIndex(x, N):
    """
    Sebastijan Mrak
    GNSS Amplitude scintillation index for the interval of the length 'N' samples
    """
    idx = np.isnan(x)
    n2 = int(N/2)
    iterate = np.arange(n2, x.size-n2)
    y = np.nan * np.copy(x)
    for i in iterate:
        if np.sum(np.isfinite(x[i-n2:i+n2])) > N/4:
            y[i] = np.sqrt(np.nanvar(x[i-n2:i+n2])) / np.nanmean(x[i-n2:i+n2])
    y[idx] = np.nan
    return y

def sigmaTEC(x, N):
    idx = np.isnan(x)
    n2 = int(N/2)
    iterate = np.arange(n2, x.size-n2)
    y = np.nan * np.copy(x)
    for i in iterate:
        if np.sum(np.isfinite(x[i-n2:i+n2])) > N/4:
            y[i] = np.nanstd(x[i-n2:i+n2])
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