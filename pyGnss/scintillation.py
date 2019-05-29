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
    roti = []    
    for i in range(len(tec)-length):
        roti.append(np.std(tec[i:i+length]))
    
    return np.array(roti)
    

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
            y[i] = np.nanstd(data[i:i+N] / np.nanmean(data[i:i+N]))
    return y

def sigmaTEC(x, N, median = False):
    idx = np.isnan(x)
    n2 = int(N/2)
    iterate = np.arange(n2, x.size-n2)
    y = np.nan * np.copy(x)
    for i in iterate:
        chunk = x[i-n2:i+n2]
#        if median:
#            x_sorted = np.sort(chunk)
#            args = np.argsort(chunk)
#            mask = np.arange(int(N/20), int(N - N/20)+1)
        if np.sum(np.isfinite(chunk)) > N/4:
            y[i] = np.nanstd(chunk)
    y[idx] = np.nan
    return y