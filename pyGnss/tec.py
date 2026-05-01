# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 11:08:09 2018

@author: smrak
"""
import numpy as np

#constnats for GPS
g1 = 1575420000
g2 = 1227600000
g5 = 1176450000
#constnats for Galieleo
e1 = 1575420000
e5 = 1176450000
e6 = 1278750000
e7 = 1207140000
e8 = 1191795000
#constants for BeiDou
c1 = 1575420000
c2 = 1561098000
c5 = 1176450000
c6 = 1268520000
c7 = 1207140000
c8 = 1191795000

c0 = 299792458.0  # Speed of light in m/s
def phase_tec(L1, L2, f1, f2):
    """
    Calculate the phase TEC from the phase measurements of two frequencies.
    
    Parameters
    ----------
    L1 : array_like
        Phase measurement at frequency f1 in cycles.
    L2 : array_like
        Phase measurement at frequency f2 in cycles.
    f1 : float
        Frequency 1 in Hz.
    f2 : float
        Frequency 2 in Hz.

    Returns
    -------
    tec : array_like
        Total Electron Content in TECU (10^16 electrons/m^2).
    """
    
    # Calculate the TEC using the formula: TEC = (delta_L * f1 * f2) / (40.3 * (f1^2 - f2^2))
    phase_tec = ((f1**2 * f2**2) / (f1**2 - f2**2)) * (c0/40.3) * \
                (L1 / f1 - L2 / f2) / pow(10, 16)
    
    return phase_tec

def range_tec(C1, C2, f1, f2):
    """
    Calculate the code TEC from the code measurements of two frequencies.
    Parameters
    ----------
    C1 : float
        Speed of light divided by frequency f1 in meters per cycle.
    C2 : float
        Speed of light divided by frequency f2 in meters per cycle.
    f1 : float
        Frequency 1 in Hz.
    f2 : float
        Frequency 2 in Hz.

    Returns
    -------
    tec : array_like
        Total Electron Content in TECU (10^16 electrons/m^2).
    """
    range_tec = ((f1**2 * f2**2) / (f1**2 - f2**2)) * (C2 - C1) / 40.3 / pow(10, 16)
    
    return range_tec