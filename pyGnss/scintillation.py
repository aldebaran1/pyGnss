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
        y[i] = np.std(data[i:i+N])
    return y
    
def AmplitudeScintillationIndex(data, N):
    """
    Sebastijan Mrak
    GNSS Amplitude scintillation index for the interval of the length 'N' samples
    """
    y = np.nan * np.zeros(data.shape[0])
    for i in range(data.shape[0] - N):
        y[i] = np.std(data[i:i+N] / np.mean(data[i:i+N]))
    return y

