#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 13:52:24 2026

@author: mraks1
"""

from pyGnss import pyGnss, tec
import xarray as xr
import georinex as gr
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np 
from scipy.signal import savgol_filter
from scipy.optimize import least_squares
from scipy import interpolate
import os, glob
from PIL import Image
import io
import xarray as xr
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
root = "/Users/mraks1/Downloads/rinex/1230/"
rnxfiles = np.random.choice(sorted(glob.glob(root + "*crx")), 1, replace=False)
for obs in rnxfiles:
# obs = '/Users/mraks1/Downloads/rinex/AJAC00FRA_S_20251660000_01D_01S_MO.crx'
# obs = '/Users/mraks1/Downloads/rinex/ALAC00ESP_S_20251660000_01D_01S_MO.crx'
    # obs = '/Users/mraks1/Downloads/rinex/1230/YIEL00AUS_S_20203650000_01D_30S_MO.crx'
    print (obs)
    # nav = '/Users/mraks1/Downloads/rinex/navtmp/GFZ0MGXRAP_20251660000_01D_05M_ORB.SP3'
    nav = '/Users/mraks1/Downloads/rinex/navtmp/GFZ0MGXRAP_20203650000_01D_05M_ORB.SP3'
    # ion = '/Users/mraks1/Downloads/rinex/ionex/ESA0OPSRAP_20251660000_01D_01H_GIM.INX'
    t0 = datetime.now()
    D = gr.load(obs, use=['G', 'E', 'C'])
    print (datetime.now()-t0)
    svlist = D.sv.values
    ts = D.interval
    if not np.isfinite(D.interval):
        ts = np.nanmedian(np.diff(D.time.values).astype('timedelta64[s]')).astype(int)
    N = gr.load(nav)
    # sat_bias = pyGnss.getSatBias(ion)
    # dcb0 = np.zeros(svlist.size)
    # if sat_bias is not None:
    #     sat_bias_keys = sat_bias.keys()
    #     for isv, sv in enumerate(svlist):
    #         if sv in sat_bias_keys:
    #             dcb0[isv] = sat_bias[sv]
    t0 = datetime.now()
    sTEC, tec_sigma, AER = pyGnss.getSTEC(D, fsp3=nav, maxgap=7, maxjump=1, el_mask=10, H=450, return_aer=True, return_tec_error=True)
    print (datetime.now()-t0)
    DCB = pyGnss.getDCBfromSTEC(sTEC, AER, el_mask=10, H=450, ts=ts)# decimate=True, ts=1)
    # DCB1 = pyGnss.getDCBfromSTEC(sTEC, AER, el_mask=5, H=450, ts=1)# decimate=True, ts=1)
    F = pyGnss.getMappingFunction(AER[:,:,1], h=450)
    
    # delta_x = np.arange(-100,100,0.5)
    # res = np.nan * np.ones((D.time.size, delta_x.size))
    # for idx, dx in enumerate(delta_x):
    #     res[:,idx] = np.nanstd(((sTEC-(dcb0+dx)) * F), axis=1)
    
    
    VTEC = (sTEC-DCB) * F
    DTEC = pyGnss.getDTECsg_from_VTEC(VTEC, N=30*60/ts, order=1)
    DTECm = np.copy(DTEC)
    DTECm[AER[:,:,1]<30] = np.nan
    mean_vtec = np.nanmean(VTEC, axis=1, keepdims=True)
    mean_vtec0 = np.nanmean(VTEC, axis=0)
    VTECm = np.copy(VTEC)
    VTECm[AER[:,:, 1]<30] = np.nan
    # ROTI = pyGnss.getROTI(sTEC=VTECm, ts=1, N=60)
    sigma_tec = pyGnss.sigmaTEC(DTEC, N=60)
    sigma_tec[AER[:,:, 1]<30] = np.nan
    
    w = np.sin(np.radians(AER[:,:,1]))**2
    
    fig = plt.figure(figsize=[10,12])
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    ax1.set_title(f"{os.path.split(obs)[1]}")
    ax1.plot(D.time, sTEC+DCB)
    ax1.set_ylabel("slant TEC")
    ax2.plot(D.time, VTECm)
    ax2.set_ylabel("vertical TEC")
    # ax3.plot(D.time, DTECm)
    ax3.plot(D.time, ((VTEC - mean_vtec)*w))#/np.sqrt(np.nanvar(VTEC, axis=1, keepdims=True)))
    # ax3.plot(D.time, mean_vtec)
    ax3.set_ylabel("dTEC")
    
    # plt.savefig(obs + ".png")
    # plt.close(fig)
    # ax4.plot(D.time, sigma_tec)
    # ax4.set_ylabel("sigma_tec")
    
    # ax4.plot(D.time, np.nanstd(VTEC*w, axis=1)**2)
    # ax4.plot(D.time, (VTEC-mean_vtec)/mean_vtec0, 'r')
    # ax4.plot(D.time, (VTEC-mean_vtec)*w, 'b')
    
    def _tec_update(delta, tec_in, ranges):
        tec_out = np.nan * np.ones(tec_in.size)
        for i, r in enumerate(ranges):
            if i == 0:
                tec_out[r[0]:r[1]] = tec_in[r[0]:r[1]]
            else:
                tec_out[r[0]:r[1]] = tec_in[r[0]:r[1]] + delta[i-1]
                
        return tec_out
    
    def _fun(x0, tec_in, ranges):
        idx = np.isfinite(tec_in)
        if np.sum(idx) != tec_in.size:
            # f = interpolate.interp1d(np.arange(tec_in.size)[idx], tec_in[idx], fill_value="extrapolate")
            # f = interpolate.interp1d(np.arange(tec_in.size)[idx], tec_in[idx], bounds_error=False, fill_value=np.nan)
            # tec_in = f(np.arange(tec_in.size))
            tec_in = np.interp(np.arange(tec_in.size), np.arange(tec_in.size)[idx], tec_in[idx])
        
        tec_out = _tec_update(x0, tec_in, ranges)
        
        res = abs(np.diff(np.diff(tec_out[np.isfinite(tec_out)])))**2
        res = res[np.isfinite(res)]
        return res
    
    frames = []
    for isv, sv in enumerate(svlist):
        if sv[0] == "G":
        # if sv == "G22":
            if int(D.version) == 3:
                L1 = D.L1C
                C1 = D.C1C
                L2 = D.L2W
                C2 = D.C2W
            else:
                L1 = D.L1
                C1 = D.C1
                L2 = D.L2
                C2 = D.P2
            f1 = g1
            f2 = g2
            # L1 = D.L1
            # C1 = D.C1
            # L2 = D.L2
            # C2 = D.P2
            # f1 = g1
            # f2 = g2
        # elif sv[0] == 'E':
        #     L1 = D.L1C
        #     C1 = D.C1C
        #     L2 = D.L5P
        #     C2 = D.C5P
        #     f1 = e1
        #     f2 = e5
        else:
            continue
        if np.sum(np.isfinite(L1.sel(sv=sv).values)) > 50:
            pTEC = tec.phase_tec(L1.sel(sv=sv).values, L2.sel(sv=sv).values, f1, f2)
            rTEC = tec.range_tec(C1.sel(sv=sv).values, C2.sel(sv=sv).values, f1, f2)
            el = AER[:,isv,1]
            idx, orbits = pyGnss.getIntervalsTEC_vec(pTEC, maxgap=10*60/ts, maxjump=None)
            # print (orbits)
            # A = pyGnss.getPhaseCorrTEC(L1.sel(sv=sv).values, L2.sel(sv=sv).values, C1.sel(sv=sv).values, C2.sel(sv=sv).values, ts=1, maxgap=1, maxjump=1.4, el=el, return_tec_err=1)
            ptec_out = np.nan * np.ones(pTEC.size)
            tec_sigma = np.nan * np.ones(len(orbits))
            for iob, orbit in enumerate(orbits):
                tmp = np.nan * np.ones(orbit[1]-orbit[0])
                _, ranges = pyGnss.getIntervalsTEC_vec(pTEC[orbit[0]: orbit[1]], maxgap=1, maxjump=1.4, min_length=10)
                
                if len(ranges) == 0:
                    continue
                elif len(ranges) > 1:
                    x0 = np.ones(len(ranges)-1)
                    sb_lsq = least_squares(_fun, x0, args=(pTEC[orbit[0]: orbit[1]], ranges), xtol=.0001)
                    x1 = sb_lsq.x
                    chi2 = 2 * sb_lsq.cost
                    tmp[:] = _tec_update(x1, pTEC[orbit[0]: orbit[1]], ranges)
                else:
                    tmp[ranges[0][0]:ranges[0][1]] = pTEC[orbit[0]:orbit[1]][ranges[0][0]:ranges[0][1]]
                
                tec_difference = tmp - rTEC[orbit[0]:orbit[1]]
                
                median_diff = np.nanmedian(tec_difference)
                
                
                w = np.sin(np.radians(AER[orbit[0]:orbit[1],isv,1]))
                offset = np.nansum(np.multiply(w, tec_difference)) / np.nansum(w)
                # print (offset, median_diff, AER[orbit[0]:orbit[1],isv,1].max())
                TECsigma = (np.nansum(np.multiply(w, tec_difference**2)) * np.nansum(w) - np.nansum(np.multiply(w,tec_difference))**2) / (np.nansum(w)**2 - np.nansum(w**2))    
                tec_sigma[iob] = np.round(np.sqrt(TECsigma),2)
                
                ptec_out[orbit[0]:orbit[1]] = tmp - median_diff
                
            
            fig = plt.figure(figsize=[8,6])
            ax1 = fig.add_subplot(211)
            ax2 = fig.add_subplot(212, sharex=ax1)
            ax1.set_title(sv)
            
            
            ax1.plot(D.time, pTEC)
            ax1.plot(D.time, ptec_out)
            
            ax2.plot(D.time, rTEC)
            ax2.plot(D.time, ptec_out, label=f"$\sigma$={tec_sigma}")
            
            ax2.legend()
    # break
            
            buf = io.BytesIO()
            plt.savefig(buf, format='png')
            buf.seek(0)
            frames.append(Image.open(buf))
            plt.close(fig)
                
    frames[0].save(f'{root} + {os.path.split(obs)[1]}.gif',
                   save_all=True,
                   append_images=frames[1:],
                   duration=500, # milliseconds
                   loop=0)         
