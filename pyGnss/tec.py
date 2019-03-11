# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 11:08:09 2018

@author: smrak
"""
import numpy as np
from typing import Union

#constnats for GPS
F1 = 1575420000
F2 = 1227600000
F5 = 1176450000
C0 = 299792458
F = {'1': F1, '2': F2, '5': F5}

def slantTEC(C1: np.ndarray = None, 
             C2: np.ndarray = None, 
             L1: np.ndarray = None, 
             L2: np.ndarray = None, 
             tec_type: str = 'range', 
             fpair: Union[list, np.ndarray] = [1, 2]):
    """
    Compute total electron content using any combination of observations. 
    Note that C1/L1 must be observations at higher frequency than C2/L2
    Specify the GPS frequency in fpair argument, [1, 2]; [1, 5]; [2, 5]
    tec in units of TECu => everything normalized by 1e16; 
    """
    K = (F[str(fpair[0])]**2 * F[str(fpair[1])]**2) / (F[fpair[0]]**2 - F[str(fpair[1])]**2)
    if tec_type == 'range':
        assert C1 is not None and C2 is not None
        tec = K / 40.3 * (C2 - C1) / 1e16
    elif tec_type == 'phase':
        assert C1 is not None and C2 is not None
        assert L1 is not None and L2 is not None
        range_tec = K / 40.3 * (C2 - C1) / 1e16
        phase_tec = K / 40.3 * C0 * (L1/F[str(fpair[0])] - L2/F[str(fpair[1])]) / 1e16
        N = np.nanmedian(range_tec - phase_tec)
        tec = phase_tec + N
    else:
        raise('TEC_TYPE must be either "phase" or "range". ')
    return tec

def getMappingFunction(elevation, H):
    """
    Return mapping function as a function of the input elevation angle and
    the height H in kilometers
    """
    
    rc1 = 6371.0 / (6371.0 + H)
    f = np.sqrt(1.0 - (np.cos(np.radians(elevation))**2 * rc1**2))
    
    return f

def getVerticalTEC(stec: np.ndarray = None,
                   F: np.ndarray = None,
                   elevation: np.ndarray = None,
                   H: Union[int, float] = None):
    """
    
    """
    if H is None and elevation is None:
        assert stec.shape[0] == F.shape[0], 'TEC and mapping function arrays must have the same length' 
    elif F is None:
        assert elevation is not None, 'Elevation abgle array must be given'
        assert H is not None, 'Ionospheric mapping function height has to be specified in kilometers'
        F = getMappingFunction(elevation, H)
    vtec = stec * F
    
    return vtec

def singleFrequencyTEC(L1: np.ndarray = None, 
                       C1: np.ndarray = None, 
                       units: str = 'm', 
                       vertical:bool = False, 
                       el: np.ndarray = np.array([]), 
                       H:int = 300,
                       frequency: int = 1):
    # L1 [cycle]: Convert C1 units to [cycles]
    if units == 'm':
        C1cycle = C1 * F1 / C0
    elif units == 'cycle':
        C1cycle = C1
    elif units == 'rad':
        C1cycle = C1 * F1 / C0 / (2 * np.pi)
    else:
        raise ('Enter appropriate units "m", "rad", or "cycle"')
    dN = np.nanmean(L1 - C1cycle)
    L1corr = L1 - dN
    TECrad = C1cycle-L1corr
    if vertical:
        assert el.shape[0] == L1.shape[0], 'observations and elevation angle arrays must be the same length'
        f = getMappingFunction(el, H)
        TECrad = TECrad * f
    
    return TECrad

def retreiveDTECfromPhase(L: np.ndarray = None, 
                          frequency: int = 1, 
                          units: str = 'cycle'):
    if units == 'cycle':
        dTEC = L * F[str(frequency)]**2 / 40.3 / 1e16
    elif units == 'rad':
        dTEC = L * C0 * F[str(frequency)] / 2 / np.pi / 40.3 /1e16
    else:
        raise ('Enter an appropriate name for units')
        
    return dTEC * -1