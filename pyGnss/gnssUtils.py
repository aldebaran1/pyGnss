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
from scipy.interpolate import CubicSpline

#constnats for GPS
f1 = 1575420000
f2 = 1227600000
f5 = 1176450000
c0 = 299792458
freq = {'1': f1, '2': f2, '5': f5}

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
def CSnorm(L1, L, cycle_slip_idx, frequency = 2, verbose = False, plot = True):
    global freq
    # Normalized values
    L1n = L1 / f1
    Ln = L / freq[str(frequency)]
    # Copy L
    Y = np.copy(L)
    # Add ix of the last entry 
    cycle_slip_idx = np.append(cycle_slip_idx, L.shape[0])
    # Zero mean strat
    N = np.nanmedian(L1n[:cycle_slip_idx[0]] - Ln[:cycle_slip_idx[0]])
    L1n -= N
    # Correct for each interval
    for i, ix in enumerate(cycle_slip_idx):
        if i == 0: continue
        ixl = cycle_slip_idx[i-1]
        if ix - ixl <= 5: continue
        Yn = Y / freq[str(frequency)]
        N = np.nanmedian(L1n[ixl:ix] - Yn[ixl:ix])
        corr = np.round(N * freq[str(frequency)])
        Y[ixl:] += corr
        if verbose: print ('{}:{}, Nnorm:{}, Corr: {}'.format(ixl,ix,N, corr))

    # Outlier interval:
    Rn = L1n - Y / freq[str(frequency)]
    stdn = np.nanstd(Rn)
    ixo = abs(Rn) > 5 * stdn
    
    if plot:
        Len = L.shape[0]
        x = np.arange(Len)
        plt.figure()
        plt.plot(Rn, 'b')
        plt.plot([0, Len], [5*stdn, 5*stdn] , 'r')
        plt.plot([0, Len], [-5*stdn, -5*stdn], 'r')
        plt.plot(x[ixo], Rn[ixo], 'xr')
    
    Y[ixo] = np.nan
    
    return Y

def phaseDetrend(y, order,polynom=False):
    """
    Sebastijan Mrak
    Raw phase data detrending using N-th polinom approximation function.
    Detrended output is input data subtracted with polinom approximation.
    Output is of the same length as input data 'y'. 
    """
    x = np.arange(y.shape[0])
    mask = np.isnan(y)
    z = np.polyfit(x[~mask], y[~mask], order)
    f = np.poly1d(z)
    polyfit = f(x)
    y_d = y-polyfit
    
    if polynom:
        return y_d, polyfit
    else:
        return y_d

def cycleSlipRepair(C, L, cycle_slip_idx, freq = 1, verbose = False, units = 'm'):
    global f1, f2, f5, c0
    if freq == 1:
        lamb = c0 / f1
    elif freq == 2:
        lamb = c0 / f2
    elif freq == 5:
        lamb = c0 / f5
    else:
        raise ('Enter the rifght frequency channel [1,2,5]')
    Y = np.copy(L)
    if freq == 1: Cd, C = phaseDetrend(C, order=12, polynom=True)
    for i, ix in enumerate(cycle_slip_idx):
        if np.isfinite(L[ix]) and np.isfinite(C[ix]):
            interval = np.arange(ix)
            finite = np.isfinite(C[interval]) & np.isfinite(Y[interval])
            if np.sum(finite) == 0: break
            last_idx = interval[finite][-1]
            if units == 'm':
                CS_m = Y[ix] - np.nanmedian(Y[last_idx-5 : last_idx])
                range_m = C[ix] - np.nanmedian(C[last_idx-5 : last_idx])
                diff = range_m - CS_m
                corr = (diff // lamb) * lamb
    #            print (corr)
                if verbose: print(ix, last_idx, CS_m, range_m, diff, corr)
            elif units == 'cycle':
                CS = Y[ix] - np.nanmedian(Y[last_idx-5 : last_idx])
                CS_L1 = C[ix] - np.nanmedian(C[last_idx-5 : last_idx])
                diff = CS_L1 - CS
                corr = np.round(diff)
                print (diff)
            else:
                break
            Y[ix:] += corr
        else:
            Y[ix] = np.nan
    return Y

def cycleSlipIdx(x, idlli):
    x_poly = phaseDetrend(x, order=12)
    mask = np.flatnonzero(np.isnan(x_poly))
    rmask = mask - 1 if mask-1 not in mask else mask - 2
    x_poly[mask] = x_poly[rmask]
    cycle_slip_idx = []
    if idlli.shape[0] > 0:
        for ix in idlli:
            if (np.sign(x_poly[ix-1]) != np.sign(x_poly[ix])) or (abs(x_poly[ix-1] - x_poly[ix]) > 1.5):
                cycle_slip_idx.append(ix)
    else:
        cycle_slip_idx = []
    return (np.array(cycle_slip_idx))



def RepairObs(D, idel, args = ['L1', 'L2', 'L5']):
    obs = {}
    for arg in args:
        argc = arg.replace('L', 'P') if ('2' in arg) else arg.replace('L', 'C')
        argi = arg + 'lli'
        try:
            C = D[argc].values[idel]
            if np.nansum(np.isfinite(C)) < 600: continue
            L = D[arg].values[idel]
            if arg == 'L1': L1 = L
            Li = D[argi].values[idel] % 2
            idlli = np.where(Li == 1)[0]
            cycle_slip_idx = cycleSlipIdx(L, idlli)
            # Convert to meters
            Lmeter = L * c0 / freq[arg[-1]]
            # Correct cycle slips
            if cycle_slip_idx.shape[0] > 0:
                print ('{} Repair'.format(argc))
                if arg == 'L1':
                    mask = np.isfinite(C)
                    x = np.arange(C.shape[0])
                    CS = CubicSpline(x[mask], C[mask])
                    C1 = CS(x)
                    L1meter = cycleSlipRepair(C1, Lmeter, cycle_slip_idx, 
                                             freq = 1, verbose = False,
                                             units = 'm')
                    L = L1meter / c0 * f1
                    L1 = L
                else:
                    L = CSnorm(L1, L, cycle_slip_idx, plot = False, verbose = False, frequency = arg[-1])
            # Convert back to cycles
            Ld, Lp = phaseDetrend(L, order = 12, polynom = True)
            x = np.arange(Ld.shape[0])
            mask= np.isnan(Lp)
            Ldp = np.copy(Ld)
            Ldp[mask] = np.interp(x[mask], x[~mask], Ld[~mask])
            Lps = hpf(Ldp)
            Lps[:15] = np.nan
        
            obs[arg] = L
            obs[argc] = C
            obs[arg+'d'] = Ldp
            obs[arg+'p'] = Lps
            obs[arg + 'lli'] = cycle_slip_idx
        except:
            pass
    return obs

def getPolynomOrder(interval, Ts, weights=[1,4,7,10]):#[1,5,9,13]
    if isinstance(interval, list):
        interval = np.array(interval)
    lst_len = (interval[1] - interval[0]) * Ts / 60 # In minutes
#    lst_len = interval.shape[0] * Ts / 60
    if lst_len < 20:
        polynom_order = weights[0]
    elif lst_len >= 20 and lst_len < 50:
        polynom_order = weights[1]
    elif lst_len >= 50 and lst_len < 100:
        polynom_order = weights[2]
    elif lst_len >= 100:
        polynom_order = weights[3]
        
#    if lst_len < 30:
#        polynom_order = weights[0]
#    elif lst_len >= 30 and lst_len < 70:
#        polynom_order = weights[1]
#    elif lst_len >= 70 and lst_len < 200:
#        polynom_order = weights[2]
#    elif lst_len >= 200:
#        polynom_order = weights[3]
    return polynom_order

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

def getPlainResidual(tec, Ts=1, maxgap=2, maxjump=1, 
                     verbose = False,typ='poly', 
                     intervals = None,
                     fc=0.005, order=3, 
                     weights=[1,4,7,10]):
    if intervals is None:
        intervals = getIntervals(tec, maxgap=maxgap, maxjump=maxjump)
    assert isinstance(intervals, list)
    tecd = np.nan * np.ones(tec.shape[0])
    for lst in intervals:
        if lst[1]-lst[0] > (60/Ts):
            polynom_order = getPolynomOrder(lst, Ts, weights=weights)
            # PolynomialOrder
#            print ('order = ', polynom_order)
            tmp = phaseDetrend(tec[lst[0]:lst[1]], order=polynom_order)
            if typ == 'hpf':
                skip = 15 if Ts > 10 else int(60/Ts)
                tmp = hpf(tmp, fs=1/Ts, fc=fc, order=order)
                tmp[:skip] = np.nan
            if np.nanmax(abs(tmp)) < 20:
                tecd[lst[0]:lst[1]] = tmp
            else:
                if verbose: print ('Something went wrong big time.')
                pass
    return tecd
