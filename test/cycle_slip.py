# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 13:59:34 2018

@author: smrak
"""
from pyGnss import pyGnss
from pyGnss import gnssUtils as gu
from pandas import Timestamp
import numpy as np
import georinex as gr
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

def mahali(folder='E:\\mahali\\',day=280,rx=9):
    # Init
    date = datetime.strptime(str(2015)+str(day), '%Y%j')
    directory = date.strftime('%Y-%m-%d')
    rxmah = {'2': 'mah2'+str(day)+'0.15o.nc',
             '3': 'mah3'+str(day)+'0.15o.nc',
             '4': 'mah4'+str(day)+'0.15o.nc',
             '5': 'mah5'+str(day)+'0.15o.nc',
             '6': 'mah6'+str(day)+'0.15o.nc',
             '7': 'mah7'+str(day)+'0.15o.nc',
             '8': 'mah8'+str(day)+'0.15o.nc',
             '9': 'mah9'+str(day)+'0.15o.nc',
             '13': 'ma13'+str(day)+'0.15o.nc',}
    folder = folder + directory + '\\'
    if isinstance(rx, int):
        nc = rxmah[str(rx)]
        fnc = folder + nc
    elif isinstance(rx, list):
        fnc = []
        for r in rx:
            nc = rxmah[str(rx)]
            fnc.append(folder + nc)
    nav = 'brdc'+str(day)+'0.15n'
    fnav = folder + nav
    
    return fnc, fnav

def TODiff(x):
    l_diff_11 = np.diff(x) # Diff 1st order
    l_diff_11 = np.hstack((np.nan, l_diff_11))
    l_diff_21 = np.diff(l_diff_11) # Diff 2nd order
    l_diff_21 = np.hstack((np.nan, l_diff_21))
    l_diff_31 = np.diff(l_diff_21) # Diff 3rd order 
    l_diff_31 = np.hstack((np.nan, l_diff_31))
    
    return l_diff_31

def getDoppler(x):
    D1 = np.diff(x) # Diff 1st order
    D1 = np.hstack((np.nan, D1))
    return D1

def detrendMetric(x, porders=[1,15]):
    value = []
    totalsum = []
    for po in porders:
        ptemp = pyGnss.phaseDetrend(x,po)
        value.append(np.sum(ptemp))
        totalsum.append(np.sum(abs(ptemp)))
    return np.array(value), np.array(totalsum)

def dataFromNC(fnc,fnav,el_mask=30,tlim=None):
    leap_seconds = gu.getLeapSeconds(fnav)
    D = gr.load(fnc, useindicators=True).sel(sv=sv)
    if tlim is not None:
        if len(tlim) == 2:
            D = D.where(np.logical_and(D.time >= np.datetime64(tlim[0]), D.time <= np.datetime64(tlim[1])), drop=True)
    obstimes64 = D.time.values
    dt = np.array([Timestamp(t).to_pydatetime() for t in obstimes64]) - \
               timedelta(seconds = leap_seconds)
    rx_xyz = D.position
    aer = pyGnss.gpsSatPosition(fnav,dt,sv=sv, rx_position=rx_xyz, coords='aer')
    idel = (aer[1] >= el_mask)
    dt = dt[idel]
    return dt, D, idel

def getIntervals(x, idx):
    ranges = []
    idx = np.hstack((0, idx, x.shape[0]))
    number_length = idx.shape[0] - 1
    if number_length > 2:
        for i in range(number_length):
            interval = [idx[i], idx[i+1]] if (i <= number_length-1) else [idx[i], idx[i+1]-1]
            interval_length = interval[1] - interval[0]
            if interval_length > 3:
                ranges.append(interval)
    else:
        ranges = [[int(0), int(x.shape[0])]]
    return ranges

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

def CSrepair(X, cycle_slip_idx):
    
    iters = 0
    if cycle_slip_idx.shape[0] > 0:
        X1 = np.copy(X)
        eps = 1e6
        Y = []
        while eps >= 1*cycle_slip_idx.shape[0] and iters < 9:
#            print (X1)
            diff_coeff = []
            iters += 1
            for i, ix in enumerate(cycle_slip_idx):
                interval = np.arange(ix)
                finite = np.isfinite(X[interval])
                last_idx = interval[finite][-1]
                last = X1[last_idx]
                Ncycle_slip = X1[ix] - last
                Ncycle_slip_norm = Ncycle_slip / (ix-last_idx)
                avg_diff = np.nanmean([X1[interval][finite][-2] - X1[interval][finite][-3],
                                       X1[interval][finite][-3] - X1[interval][finite][-4],
                                       X1[interval][finite][-4] - X1[interval][finite][-5],
                                       X1[interval][finite][-5] - X1[interval][finite][-6]])
                corr = round(Ncycle_slip_norm - avg_diff)
                diff_coeff.append(abs(corr))
                if ix == cycle_slip_idx[-1]:
                    X1[ix:] = X1[ix:] - corr
                else:
                    X1[ix : cycle_slip_idx[i+1]] = X1[ix : cycle_slip_idx[i+1]] - corr
            print(diff_coeff)
            eps = np.sum(diff_coeff)
            Y.append(phaseDetrend(X1, order=12))
    else:
        iters+=1
        Y = [phaseDetrend(X, order=12)]
    
    return np.array(Y), iters

def dydx(x):
    DYDX = np.nan * np.copy(x)
    mask = np.isfinite(x)
    x1 = np.arange(x.shape[0])
    X = x1[mask]
    Y = x[mask]
    dy = np.hstack((np.nan, np.diff(Y) / np.diff(X)))
    DYDX[mask] = dy
    return DYDX


def CSrepair2(X, cycle_slip_idx):
    
    iters = 0
    if cycle_slip_idx.shape[0] > 0:
        X1 = np.copy(X)
        eps = np.ones(cycle_slip_idx.shape[0]) * 1e3
        Y = []
        while abs(eps).max() >= 1 and iters < 15:
            iters += 1
            for i, ix in enumerate(cycle_slip_idx):
                interval = np.arange(ix)
                finite = np.isfinite(X[interval])
                last_idx = interval[finite][-1]
                #derivative
                dy = dydx(X1)
#                plt.plot(dy)
                last = dy[last_idx]
                Ncycle_slip = dy[ix] - last
                # Correction step
                dx = (ix-last_idx)
                Ncycle_slip_norm = Ncycle_slip / dx
                avg_diff = np.nanmean([dy[interval][finite][-2] - dy[interval][finite][-3],
                                       dy[interval][finite][-3] - dy[interval][finite][-4],
                                       dy[interval][finite][-4] - dy[interval][finite][-5],
                                       dy[interval][finite][-5] - dy[interval][finite][-6]])
                Nstep_nominal = avg_diff * dx
                
                if abs(Ncycle_slip_norm) > 1:
                    
                    corr = np.ceil(Ncycle_slip_norm - avg_diff) if (Ncycle_slip_norm - avg_diff) > 2 else np.floor(Ncycle_slip_norm - avg_diff)
                    eps[i] = corr
                    X1[ix:] = X1[ix:] - corr
                    print (ix, avg_diff, Ncycle_slip, corr)
                elif abs(Nstep_nominal - Ncycle_slip) > 1 and abs(Ncycle_slip_norm) <= 1:
                    corr = np.round(Ncycle_slip - Nstep_nominal)
                    eps[i] = corr
                    X1[ix:] = X1[ix:] - corr
                    print (ix, Nstep_nominal, Ncycle_slip, Nstep_nominal - Ncycle_slip, corr)
#                    print(corr)
#                    print (ix, corr)
#                    if ix == cycle_slip_idx[-1]:
#                        X1[ix:] = X1[ix:] - corr
#                    else:
#                        X1[ix : cycle_slip_idx[i+1]] = X1[ix : cycle_slip_idx[i+1]] - corr
                else:
                    eps[i] = 0
#                print (Ncycle_slip_norm, avg_diff)
#            print(eps)
#            eps = np.sum(diff_coeff)
            Y.append(phaseDetrend(X1, order=12))
    else:
        iters+=1
        Y = [phaseDetrend(X, order=12)]
        X1 = X
    
    return X1, np.array(Y), iters
###############################################################################
day = 280
year = 2015
tlim = [datetime.strptime("{} {} {}:{}".format(year,day,6,5), "%Y %j %H:%M"), 
        datetime.strptime("{} {} {}:{}".format(year,day,11,5), "%Y %j %H:%M")]
folder = 'C:\\Users\\smrak\\Google Drive\\BU\\projects\\Themahali\\data\\rinex\\'
args = ['C1', 'L1', 'L2', 'S1', 'vTEC']
arg = 'L1'
sv = 'G23'
rx = 8
fnc, fnav = mahali(folder = folder, rx = rx, day = day)

dt, D, idel = pyGnss.dataFromNC(fnc,fnav,sv=sv,tlim=tlim)

# L1
C1 = D['C1'].values[idel]
L1 = D['L1'].values[idel]
L1i = D['L1lli'].values[idel]
L1polyfit = pyGnss.phaseDetrend(L1, order=12)
D1 = getDoppler(L1)
#plt.plot(L1polyfit)


# L2
C2 = D['P2'].values[idel]
L2 = D['L2'].values[idel]
L2i = D['L2lli'].values[idel] % 2 # Antispoofing = 4, LLi = 1,5,7 ie, modulo 2 remeinder
D2 = getDoppler(L2)
lli_sum = np.nansum(L2i)
print (np.nansum(L2i))
idlli = np.where(L2i == 1)[0]
lli_time = dt[idlli]

L2polyfit = phaseDetrend(L2, order=12)
cycle_slip_idx = cycleSlipIdx(L2, idlli)
ranges = getIntervals(L2polyfit, cycle_slip_idx)
# derivative
DYDX_1 = dydx(L2)

#dydx = np.nan * np.copy(L2)
#for r in ranges:
#    interval = np.arange(r[0], r[1])
#    finite = np.isfinite(L2[r[0] : r[1]])
#    dydx = np.diff(L2[interval][finite]) / np.diff(interval[finite])
#plt.plot(dydx)

L2, L2corrpolyfit, iters = CSrepair2(L2, cycle_slip_idx)
DYDX_2 = dydx(L2)
#iters = 0
#if cycle_slip_idx.shape[0] > 0:
#    L2corr = np.copy(L2)
#    eps = 1e6
#    L2corrpolyfit = []
#    while eps >= 1*cycle_slip_idx.shape[0] and iters < 9:
#        diff_coeff = []
#        iters += 1
#        for i, ix in enumerate(cycle_slip_idx):
#            interval = np.arange(ix)
#            finite = np.isfinite(L2[interval])
#            last_idx = interval[finite][-1]
#            last = L2corr[last_idx]
#            Ncycle_slip = L2corr[ix] - last
#            Ncycle_slip_norm = Ncycle_slip / (ix-last_idx)
#            avg_diff = np.nanmean([L2[interval][finite][-2] - L2[interval][finite][-3],
#                               L2[interval][finite][-3] - L2[interval][finite][-4],
#                               L2[interval][finite][-4] - L2[interval][finite][-5],
#                               L2[interval][finite][-5] - L2[interval][finite][-6]])
#            corr = round(Ncycle_slip_norm - avg_diff)
#            diff_coeff.append(abs(corr))
#    #        print (corr)
#            if ix == cycle_slip_idx[-1]:
#                L2corr[ix:] = L2corr[ix:] - corr
#            else:
#                L2corr[ix : cycle_slip_idx[i+1]] = L2corr[ix : cycle_slip_idx[i+1]] - corr
#        print(diff_coeff)
#        eps = np.sum(diff_coeff)
#        L2corrpolyfit.append(phaseDetrend(L2corr, order=12))
#else:
#    iters+=1
#    L2corrpolyfit.append(phaseDetrend(L2corr, order=12))

fig = plt.figure(figsize=(10,6))
plt.plot(dt, L1polyfit, 'r', lw=2)
plt.plot(dt, L2corrpolyfit[-1,:], 'b', lw=2)
plt.scatter(dt[cycle_slip_idx], L2corrpolyfit[-1,cycle_slip_idx],
            marker='x',color='k', s=200)
date_formatter = DateFormatter('%H:%M')
ax = fig.add_subplot(111)
#for i in range(iters):
#    ax.plot(dt, L2corrpolyfit[i,:], label=str(i+1))
plt.legend()
#plt.plot(dt, L1polyfit, 'r', lw=2)
ax.xaxis.set_major_formatter(date_formatter)