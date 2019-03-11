#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 14:49:04 2017

@author: Sebastijan Mrak <smrak@gmail.com>
"""
import numpy as np
from scipy import signal
import datetime
import matplotlib.pyplot as plt

# %% Filtering
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
#    mid = ((highcut - lowcut)/3)*2 / nyq
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
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Group delay [samples]')
        plt.margins(0, 0.1)
        plt.grid(which='both', axis='both')
    return y_filt , gd[idx]

def hpf(y, fc=0.1, order=5, fs=1,plot=False, group_delay=False, verbatim=False):
    """
    Sebastijan Mrak
    Filter the input data 'y' with desired HP filter.  
    """
    b, a = butter_hpf(fc, fs, order)
    y_filt = signal.lfilter(b, a, y)
    w, h = signal.freqz(b, a, worN=1000)
    gd = -np.diff(np.unwrap(np.angle(h)))/np.diff(w)
    if verbatim: print ('Group delay of the filter is '+ str(gd[-1])+' samples.')
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
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Group delay [samples]')
        plt.margins(0, 0.1)
        plt.grid(which='both', axis='both')
    if group_delay:
        return y_filt, gd[-1]
    else:
        return y_filt

def lpf(y, fc=0.1, order=5, fs=1, plot=False, group_delay=False, verbatim=False):
    b, a = butter_lpf(fc, fs, order)
    y_filt = signal.lfilter(b, a, y)
    w, h = signal.freqz(b, a, worN=1000)
    gd = -np.diff(np.unwrap(np.angle(h)))/np.diff(w)
    if verbatim: print ('Group delay of the filter is '+ str(gd[1])+' samples.')
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
    
# %% 
def datetime2posix(dtime):
    """
    Convert an input list of datetime format timestamp to posix timestamp
    """
    return [i.replace(tzinfo=datetime.timezone.utc).timestamp() for i in dtime]
# %%
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

# %%
def hilbertTransform(x, fs=1):
    """
    Return Analytic signal (envelope), phase and frequency from the Hilber Transform
    """
    from scipy.signal import hilbert
    if np.nansum(np.isnan(x)) > 0:
        y = np.nan * np.copy(x)
        mask = np.isfinite(x)
        x = x[mask]
    analytic_signal = hilbert(x)
    amplitude_envelope = np.abs(analytic_signal)
    instantaneous_phase = np.unwrap(np.angle(analytic_signal))
    instantaneous_frequency = (np.diff(instantaneous_phase) /
                               (2.0*np.pi) * fs)
    if 'y' in locals():
        amp = np.copy(y)
        amp[mask] = amplitude_envelope
        phase = np.copy(y)
        phase[mask] = instantaneous_phase
        frequency = np.copy(y)
        mask[0] = False
        frequency[mask] = instantaneous_frequency
        return amp, phase, frequency
    else:
        return amplitude_envelope, instantaneous_phase, instantaneous_frequency
    
# %%
        