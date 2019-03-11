# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 13:59:34 2018

@author: smrak
"""
import sys
import xarray
from pyGnss import pyGnss
from pyGnss import gnssUtils as gu
from pandas import Timestamp
import dascutils as du
import numpy as np
import georinex as gr
from datetime import datetime, timedelta
from scipy.interpolate import CubicSpline
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
             '13': 'ma13'+str(day)+'0.15o.nc',
             'clgo': 'clgo'+str(day)+'0.15o.nc',
             'ab33': 'ab33'+str(day)+'0.15o.nc'}
    folder = folder + directory + '\\'
    if isinstance(rx, (int, str)):
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
    idx = np.hstack((0, idx, x.shape[0]-1))
    number_length = idx.shape[0] - 1
    if number_length > 2:
        for i in range(number_length):
            interval = [idx[i], idx[i+1]] if (i <= number_length-1) else [idx[i], idx[i+1]-1]
            interval_length = interval[1] - interval[0]
            if interval_length > 3:
                ranges.append(interval)
    else:
        ranges = [[int(0), int(x.shape[0])-1]]
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
    X1 = np.copy(X)
    iters = 0
    if cycle_slip_idx.shape[0] > 0:
        eps = np.ones(cycle_slip_idx.shape[0]) * 1e3
        Y = []
        while abs(eps).max() >= 1 and iters < 1:
            iters += 1
            for i, ix in enumerate(cycle_slip_idx):
                interval = np.arange(ix)
                finite = np.isfinite(X1[interval])
                last_idx = interval[finite][-1]
                #derivative
                dy = dydx(X1)
                last = dy[last_idx]
                Ncycle_slip = dy[ix] - last
                # Correction step
                dx = (ix-last_idx)
                Ncycle_slip_norm = Ncycle_slip / dx
                avg_diff = np.nanmean([dy[interval][finite][-2] - dy[interval][finite][-3],
                                       dy[interval][finite][-3] - dy[interval][finite][-4],
                                       dy[interval][finite][-4] - dy[interval][finite][-5],
                                       dy[interval][finite][-5] - dy[interval][finite][-6]])
                avg_diff_2 = np.nanmean([dy[ix+1] - dy[ix],
                                         dy[ix+2] - dy[ix+1],
                                         dy[ix+3] - dy[ix+2],
                                         dy[ix+4] - dy[ix+3]])
                tracking_step = np.mean([avg_diff, avg_diff_2]) * dx
                
                if abs(Ncycle_slip_norm) > 1:
                    corr = np.ceil(Ncycle_slip_norm - avg_diff) if (Ncycle_slip_norm - avg_diff) > 2 else np.floor(Ncycle_slip_norm - avg_diff)
                    eps[i] = corr
                    X1[ix:] = X1[ix:] - corr
                    print (ix, avg_diff, Ncycle_slip, corr)
                elif abs(tracking_step - Ncycle_slip) > 1 and abs(Ncycle_slip_norm) <= 1:
                    corr = np.round(tracking_step - Ncycle_slip)
                    eps[i] = corr
                    X1[ix:] += (corr*dx)
                    print (ix, tracking_step, dx, Ncycle_slip, corr)
#                    print(corr)
#                    print (ix, corr)
#                    if ix == cycle_slip_idx[-1]:
#                        X1[ix:] = X1[ix:] - corr
#                    else:
#                        X1[ix : cycle_slip_idx[i+1]] = X1[ix : cycle_slip_idx[i+1]] - corr
                else:
                    eps[i] = 0
                plt.plot(phaseDetrend(X1, order=12)[2500:5000])
                
#                print (Ncycle_slip_norm, avg_diff)
#            print(eps)
#            eps = np.sum(diff_coeff)
            Y.append(phaseDetrend(X1, order=12))
    else:
        iters+=1
        Y = [phaseDetrend(X, order=12)]
        X1 = X
    
    return X1, np.array(Y), iters

def CSrepair3(X, cycle_slip_idx):
    X1 = np.copy(X)
    iters = 0
    if cycle_slip_idx.shape[0] > 0:
        eps = np.ones(cycle_slip_idx.shape[0]) * 1e3
        Y = []
        while abs(eps).max() >= 1 and iters < 5:
            iters += 1
            for i, ix in enumerate(cycle_slip_idx):
                interval = np.arange(ix)
                finite = np.isfinite(X1[interval])
                last_idx = interval[finite][-1]
                #Plot
#                plt.plot(phaseDetrend(X1, order=12))
                #derivative
                y = phaseDetrend(X1, order=12)
                last = y[last_idx]
                dx = ix - last_idx
                CS = y[ix] - last
#                print (ix, CS)
                # Correction step
#                Ncycle_slip_norm = Ncycle_slip / dx
                if abs(CS) > 1:
                    corr = np.round(CS / dx)
                    eps[i] = corr
                    X1[ix:] = X1[ix:] - corr
#                    print (ix, last, y[ix], CS, corr)

#                if abs(Ncycle_slip_norm) > 1:
#                    corr = np.ceil(Ncycle_slip_norm - avg_diff) if (Ncycle_slip_norm - avg_diff) > 2 else np.floor(Ncycle_slip_norm - avg_diff)
#                    eps[i] = corr
#                    X1[ix:] = X1[ix:] - corr
#                    print (ix, avg_diff, Ncycle_slip, corr)
#                elif abs(tracking_step - Ncycle_slip) > 1 and abs(Ncycle_slip_norm) <= 1:
#                    corr = np.round(tracking_step - Ncycle_slip)
#                    eps[i] = corr
#                    X1[ix:] += (corr*dx)
#                    print (ix, tracking_step, dx, Ncycle_slip, corr)
#                    print(corr)
#                    print (ix, corr)
#                    if ix == cycle_slip_idx[-1]:
#                        X1[ix:] = X1[ix:] - corr
#                    else:
#                        X1[ix : cycle_slip_idx[i+1]] = X1[ix : cycle_slip_idx[i+1]] - corr
                else:
                    eps[i] = 0
#                plt.plot(phaseDetrend(X1, order=12))
#                print (iters)
#                print (Ncycle_slip_norm, avg_diff)
#            print(eps)
#            eps = np.sum(diff_coeff)
            Y.append(phaseDetrend(X1, order=12))
    else:
        iters+=1
        Y = [phaseDetrend(X, order=12)]
        X1 = X
    
    return X1, np.array(Y), iters

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
    
def phaseTEC(C1,C2,L1,L2, cycle_slip_idx):
    global f1, f2
    intervals = getIntervals(L2, cycle_slip_idx)
    tec = np.nan * np.copy(L1)
    rangetec = np.nan * np.copy(L1)
    for r in intervals:
        interval = np.arange(r[0], r[1]+1)
#        print (interval)
        if r[1] - r[0] > 10:
            rtec = (f1**2 * f2**2) / (f1**2 - f2**2) / 40.3 * (C2[interval] - C1[interval]) / 1e16
            ptec = (f1**2 * f2**2) / (f1**2 - f2**2) / 40.3 * (L1[interval] - L2[interval]) / 1e16
            corr = np.nanmedian(rtec) - np.nanmedian(ptec)
#            print (corr)
            tec[interval] = ptec + corr
            rangetec[interval] = rtec
    return tec, rangetec

def slantTEC(C1, C2, L1, L2, frequency = 2):
    global freq, c0
    F = (freq['1']**2 * freq[str(frequency)]**2) / (freq['1']**2 - freq[str(frequency)]**2)
    rangetec = F / 40.3 * (C2 - C1) / 1e16
    phasetec = F / 40.3 * c0 * (L1/freq['1'] - L2/freq[str(frequency)]) / 1e16
    N = np.nanmedian(rangetec - phasetec)
    return phasetec + N
    
def getMappingFunction(elevation, H):
    rc1 = 6371.0 / (6371.0 + 200)
    f = np.sqrt(1 - (np.cos(np.radians(elevation))**2 * rc1**2))
    return f

def _processTEC(obs, sv, frequency = 2):
    stec = slantTEC(obs['C1'], obs['P2'], obs['L1'], obs['L2'], frequency = frequency)
    stec += satbias[sv]
    elevation = D.el.values
    f = getMappingFunction(elevation, H = 200)
    vtec = stec * f
    tecd = phaseDetrend(vtec, order=12)
    x = np.arange(tecd.shape[0])
    mask= np.isnan(tecd)
    tecdp = np.copy(tecd)
    tecdp[mask] = np.interp(x[mask], x[~mask], tecd[~mask])
    tecps = gu.hpf(tecdp)
    tecps[:15] = np.nan
    
    return vtec, tecd, tecps

def RepairObs(D, args = ['L1', 'L2', 'L5']):
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
            Lps = gu.hpf(Ldp)
            Lps[:15] = np.nan
        
            obs[arg] = L
            obs[argc] = C
            obs[arg+'d'] = Ldp
            obs[arg+'p'] = Lps
            obs[arg + 'lli'] = cycle_slip_idx
        except:
            pass
    return obs
f1 = 1575420000
f2 = 1227600000
f5 = 1176450000
c0 = 299792458
freq = {'1': f1, '2': f2, '5': f5}
###############################################################################
day = 280
year = 2015
tlim = [datetime.strptime("{} {} {}:{}".format(year,day,6,5), "%Y %j %H:%M"), 
        datetime.strptime("{} {} {}:{}".format(year,day,15,5), "%Y %j %H:%M")]
folder = 'C:\\Users\\smrak\\Google Drive\\BU\\projects\\Themahali\\data\\rinex\\'
savefolder = 'C:\\Users\\smrak\\Google Drive\\BU\\projects\\Themahali\\plots\\phase_detrend\\'+str(day)+'\\'

satbias = pyGnss.getSatBias(folder+datetime.strftime(tlim[0], '%Y-%m-%d')+'\\')

dascfn = 'C:\\Users\\smrak\\Google Drive\\BU\\projects\\Themahali\\data\\dasc\\' + datetime.strftime(tlim[0], '%Y-%m-%d')+'\\'
# DASC
dasc_data = xarray.open_dataset(dascfn+'interp.nc', group='DASC', autoclose = True)


args = ['L1', 'L2', 'L5']
rx = [2, 3, 4, 5, 6, 7, 8, 9, 13, 'clgo', 'ab33']
rx = ['clgo', 'ab33']
el_mask = 30
save = True
for rx in rx:
    try:
        fnc, fnav = mahali(folder = folder, rx = rx, day = day)
        svlist = gr.load(fnc).sv.values
        navdata = gr.load(fnav)
        navdatatime = navdata.time.values

        for sv in svlist:
#            sv = 'G09'
            try: 
                dt, D, idel = pyGnss.dataFromNC(fnc,fnav,sv=sv,tlim=tlim,el_mask=el_mask, satpos=True, ipp=True, ipp_alt = 110)
            except: continue
            if np.nansum(D['C1'].values[idel]) < 500: 
                continue
            try: obs = RepairObs(D, args)
            except Exception as e:
                print (e)
                continue
            #TEC
            try: tec, tecd, tecps = _processTEC(obs, sv, frequency = 2)
            except Exception as e:
                print (e)
                continue
            #LOS Intercept Brightness
            try:
                ipp_lon = D.ipp_lon.values
                ipp_lat = D.ipp_lat.values
                dasc_time, intercept_brightenss = du.getPixelBrightness(dasc_data, obstimes=dt, obs_lat=ipp_lat, obs_lon=ipp_lon)
                idT = np.array([abs(dasc_time - it).argmin() for it in dt])
                dasc_rescale = intercept_brightenss[idT]
            except Exception as e:
                print (e)
                continue
            # FFT
            NFFT = 2**11
            Tp = 1
            Lomega = np.fft.fft(obs['L1d'][np.isfinite(obs['L1d'])], n=NFFT)
            N = Lomega.shape[0]
            SLomega = abs(Lomega[:int(N/2)])**2
            SLdb = np.log10(SLomega)
            SLdbn = SLdb - np.nanmax(SLdb)
            Lfreq = np.fft.fftfreq(N, d=Tp)[:int(N/2)]
            TEComega = np.fft.fft(tecd[np.isfinite(tecd)], n=NFFT)
            STEComega = abs(TEComega[:int(N/2)])**2
            STECdb = np.log10(STEComega)
            STECdbn = STECdb - np.nanmax(STECdb)
            
###############################################################################
            #Figure
            fig = plt.figure(figsize=(18,12))
            date_formatter = DateFormatter('%H:%M')
            try:
                ax = fig.add_subplot(431)
                ax.plot(dt, obs['L1d'], 'b', lw=2, label='L1')
                ax.plot(dt, obs['L2d'], 'r', lw=2, label='L2')
                if obs['L2lli'].shape[0] > 0:
                    plt.scatter(dt[obs['L2lli']], obs['L2d'][obs['L2lli']],
                                marker='x',color='g', s=200)
                ax.set_ylabel('$\delta \Phi [m]$')
                ax.grid(axis='y')
                ax.xaxis.set_major_formatter(date_formatter)
            except: pass
            # TEC
            try:
                ax2 = fig.add_subplot(432)
                ax2.grid(axis='y')
                ax2.set_ylabel('vTEC [TECu]')
                ax2.plot(dt, tec, 'k')
                ax22 = ax2.twinx()
                ax22.plot(dt, tecd, 'b')
                ax2.xaxis.set_major_formatter(date_formatter)
            except: pass
            # Brightness
            try:
                ax3 = fig.add_subplot(433)
                ax3.grid(axis='y')
                ax3.plot(dasc_time, intercept_brightenss, '.g')
                ax3.plot(dasc_time, intercept_brightenss, 'g')
                ax3.xaxis.set_major_formatter(date_formatter)
                ax3.set_ylabel('brightness')
            except: pass
            #Scint phase
            try:
                ax4 = fig.add_subplot(434)
                ax4.plot(dt, obs['L2p'], '.r', label='L2')
                ax4.plot(dt, obs['L1p'], '.b', label='L1')
                ax4.xaxis.set_major_formatter(date_formatter)
                ax4.grid(axis='y')
                ax4.set_ylabel('$hpf(\delta \Phi) [m]$')
                ax4.legend()
            except: pass
            # ROT
            try:
                ax5 = fig.add_subplot(435)
                ax5.plot(dt, tecps, '.k')
                ax5.grid(axis='y')
                ax5.set_ylabel('$hpf(\delta$ TEC) [TECu]')
                ax5.xaxis.set_major_formatter(date_formatter)
            except: pass
            # Corr TEC, brightnes
            try:
                ax6 = fig.add_subplot(4,3,6)
                ax6.plot(dasc_rescale, abs(tecd), '.b', label='558_tec')
                Xm = np.nanmax(dasc_rescale)
                Ym = np.nanmax(abs(tecd))
                ax6.text(Xm - 0.2*Xm, Ym - 0.3*Ym,
                         s='$C_ij$ = {0:.2f}'.format(np.corrcoef(dasc_rescale[np.isfinite(tec)], tecd[np.isfinite(tecd)])[0,1]))
                ax6.set_xlabel('brightenss')
                ax6.set_ylabel('$\delta$ TEC')
                ax6.legend()
            except: pass
            # SNR
            try:
                ax7 = fig.add_subplot(437)
                ax7.plot(dt, D['S1'].values[idel], 'b', label = 'S1')
                ax7.plot(dt, D['S2'].values[idel], '--r', label = 'S2')
                if 'L5' in obs: 
                    ax.plot(dt, obs['L5d'], 'g', lw=2, label='L5')
                    ax7.plot(dt, D['S5'].values[idel], '-.g', label = 'S5')
                ax7.set_ylabel('SNR [dB]')
                ax.legend()
                ax7.legend()
                ax7.xaxis.set_major_formatter(date_formatter)
            except: pass
            # PSD
            try:
                ax8 = fig.add_subplot(438)
                ax8.plot(Lfreq, SLdbn, '.b', label='L1')
                ax8.plot(Lfreq, STECdbn, '.k', label='TEC')
                ax8.legend()
                ax8.set_ylabel('PSD [dB/Hz]')
                ax8.set_xlabel('Freq [Hz]')
            except: pass
            # Corr TEC, brightnes
            try:
                ax9 = fig.add_subplot(4,3,9)
                ax9.plot(dasc_rescale, abs(obs['L1d']), '.b', label='558_$\delta TEC$')
                Xm = np.nanmax(dasc_rescale)
                Ym = np.nanmax(abs(obs['L1d']))
                ax9.text(Xm - 0.2*Xm, Ym - 0.3*Ym,
                         s='$C_ij$ = {0:.2f}'.format(np.corrcoef(dasc_rescale[np.isfinite(obs['L1d'])], obs['L1d'][np.isfinite(obs['L1d'])])[0,1]))
                ax9.set_xlabel('brightenss')
                ax9.set_ylabel('$\delta \Phi$')
                ax9.legend()
            except: pass
            # Corr
            try:
                ax10 = fig.add_subplot(4,3,10)
                ax10.plot(abs(obs['L1p']), abs(tecps), '.b', label='hpf, $f_c=$0.1')
                Xm = np.nanmax(abs(obs['L1p']))
                Ym = np.nanmax(abs(tecps))
                ax10.text(Xm - 0.2*Xm, Ym - 0.3*Ym,
                         s='$C_ij$ = {0:.2f}'.format(np.corrcoef(-obs['L1p'][np.isfinite(obs['L1p'])], tecps[np.isfinite(obs['L1p'])])[0,1]))
                ax10.set_xlabel('$\delta$ L1')
                ax10.set_ylabel('$\delta$ TEC')
                ax10.legend()
            except: pass
            # Corr
            try:
                ax11 = fig.add_subplot(4, 3, 11)
                ax11.plot(abs(obs['L1d']), abs(tecd), '.b')
                Xm = np.nanmax(abs(obs['L1d']))
                Ym = np.nanmax(abs(tecd))
                ax11.text(Xm - 0.2 * Xm, Ym - 0.3 * Ym,
                         s='$C_ij$ = {0:.2f}'.format(np.corrcoef(-obs['L1d'][np.isfinite(obs['L1d'])], tecd[np.isfinite(obs['L1d'])])[0,1]))
                ax11.set_xlabel('$\delta$ L1')
                ax11.set_ylabel('$\delta$ TEC')
            except: pass
            # SNR, brightnes
            try:
                ax12 = fig.add_subplot(4,3,12)
                ax12.plot(dasc_rescale, D['S1'].values[idel], '.b', label='558_S1')
                Xm = np.nanmax(dasc_rescale)
                Ym = np.nanmax(D['S1'].values[idel])
                ax12.text(Xm - 0.2 * Xm, Ym - 0.3 * Ym,
                         s='$C_ij$ = {0:.2f}'.format(np.corrcoef(dasc_rescale[np.isfinite(D['S1'].values[idel])], D['S1'].values[idel][np.isfinite(D['S1'].values[idel])])[0,1]))
                ax12.set_xlabel('brightenss')
                ax12.set_ylabel('SNR S1$')
                ax12.legend()
            except: pass
            fig.suptitle('Receiver: {}; SV: {}; Day: {}'.format(rx, sv, day))
            # Save
            try:
                if save:
                    figtitle = savefolder + '{}_mah{}_{}.png'.format(day, rx, sv)
                    plt.savefig(figtitle, dpi=300)
                    plt.close(fig)
                else:
                    plt.show()
            except Exception as e: 
                exc_type, exc_obj, exc_tb = sys.exc_info()
                print (e + 'Lin: ' + exc_tb.tb_lineno)
#            break
    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        print ('Lin: ' + str(exc_tb.tb_lineno))
#    break
