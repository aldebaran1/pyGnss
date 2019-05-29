#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 10:34:49 2019

@author: smrak
"""

import numpy as np
from scipy import signal
from datetime import datetime, timedelta
from typing import Union
from scipy import ndimage
import matplotlib.pyplot as plt

def getNeighbours(image,i,j, N=1, direction='both', mask=False):
    """
    Return an array of <=9 neighbour pixel of an image with a center at (i,j)
    """
    nbg = []
    nd = N
    nu = N + 1
    if direction in['x',  'lon']:
        irange = np.arange(i - nd, i + nu)
        jrange = [j]
    elif direction in ['y', 'lat']:
        irange = [i]
        jrange = np.arange(j - nd, j + nu)
    else:
        irange = np.arange(i - nd, i + nu)
        jrange = np.arange(j - nd, j + nu)
    for k in irange:
        for l in jrange:
            try:
                nbg.append(image[k,l])
            except Exception as e:
                pass
    if mask:
        return np.array(nbg), irange, jrange
    else:
        return np.array(nbg)
    
def fillPixels(im, N=1):
    """
    Fill in the dead pixels. If a dead pixel has a least 4 finite neighbour
    pixel, than replace the center pixel with a mean valuse of the neighbours
    """
    im0 = np.copy(im)
    for n in range(N):
        skip = int(np.floor((2+n)/2))
        iterate_j_f = np.arange(0,im.shape[1],skip) if (n%2 == 0) else np.arange(im.shape[1]-1,-1,-skip)
        iterate_j_b = np.arange(0,im.shape[1],skip) if (n%2 == 1) else np.arange(im.shape[1]-1,-1,-skip)
        iterate_i_f = np.arange(0, int(im.shape[0] * 0.6), skip) if (n % 2 == 0) else np.arange(im.shape[0]-1, int(im.shape[1] * 0.4), -skip)
        iterate_i_b = np.arange(0, int(im.shape[1] * 0.6), skip) if (n % 2 == 1) else np.arange(im.shape[1]-1, int(im.shape[1] * 0.4), -skip)
        for i in iterate_i_f:
            for j in iterate_j_f:
                # Check if the pixel is dead, i.e. empty
                if np.isnan(im[i,j]):
                    # Get its neighbours as a np array
                    nbg = getNeighbours(im0,i,j,N=(3+n))
                    # If there are at leas 4 neighbours, replace the value with a mean
                    if sum(np.isfinite(nbg)) >= 4:
                        im[i,j] = np.nanmean(nbg)
        for i in iterate_i_b:
            for j in iterate_j_b:
                # Check if the pixel is dead, i.e. empty
                if np.isnan(im[i,j]):
                    # Get its neighbours as a np array
                    nbg = getNeighbours(im0,i,j,N=(3+n))
                    # If there are at leas 4 neighbours, replace the value with a mean
                    if sum(np.isfinite(nbg)) >= 4:
                        im[i,j] = np.nanmean(nbg)
    return im

def keogram(time: np.ndarray=None,
            xgrid: np.ndarray=None,
            ygrid: np.ndarray=None,
            im: np.ndarray=None,
            Xt: Union[int,float,list,np.ndarray] = None,
            Yt: Union[int,float,list,np.ndarray] = None,
            tlim: Union[list, np.ndarray] = None,
            line: Union[int, bool] = False,
            Xn: int = 0,
            Yn: int = 0,
            fillPixel: int = None,
            plotmask: bool = False,
            imskip: int = None,
            filter:str = None):
    if imskip is None: imskip = 1
    # Make a grid of input coordinates
    assert len(xgrid.shape) == len(ygrid.shape)
    if len(xgrid.shape) == 1:
        X,Y = np.meshgrid(ygrid,xgrid)
    else:
        X = xgrid
        Y = ygrid
    # Filter by time
    if tlim is not None:
        if isinstance(tlim, list):
            assert len(tlim) == 2
        elif isinstance(tlim, datetime):
            tlim = [tlim, tlim]
        else:
            raise ('Something went wrong in treq indexing')
        idT = (time >= tlim[0]) & (time <= tlim[1])
        time = time[idT]
        im = im[idT]
        im = im[::imskip]
    # Fill Pixel:
    if fillPixel is not None:
        assert isinstance(fillPixel, int)
        im[0] = fillPixels(np.squeeze(im), N=fillPixel)
    # Time dimension
    if isinstance(time, (list, np.ndarray)):
        tlen = int(time.shape[0] / imskip)
    else:
        tlen = 1
    # Filter in space
    if isinstance(Xt, (int, float)):
        idX = (xgrid >= Xt) & (xgrid <= Xt)
        if np.sum(idX) == 0: 
            idX = np.zeros((xgrid.size), dtype=bool)
            ix = abs(xgrid - Xt).argmin()
            idX[ix] = True
    elif isinstance(Xt, (list, np.ndarray)):
        idX = (xgrid >= Xt.min()) & (xgrid <= Xt.max())
    if isinstance(Yt, (int, float, np.int64)):
        idY = (ygrid >= Yt) & (ygrid <= Yt)
        if np.sum(idY) == 0: 
            idY = np.zeros((ygrid.size), dtype=bool)
            iy = abs(ygrid - Yt).argmin()
            idY[iy] = True
    elif isinstance(Yt, (list, np.ndarray)):
        idY = (ygrid >= Yt.min()) & (ygrid <= Yt.max())
    # Filter in space if doing for an arbitrary line
    if line == True:
        mask = np.full((xgrid.shape[0], ygrid.shape[0]), False, dtype = bool)
        for i in range(Xt.shape[0]):
            ix = abs(xgrid - Xt[i]).argmin()
            iy = abs(ygrid - Yt[i]).argmin()
            if not mask[ix,iy]: 
                mask[ix,iy] = True
    else:
        mask = np.meshgrid(idY,idX)[0] & np.meshgrid(idY,idX)[1]
    if np.sum(idX) > np.sum(idY):
        Yaxis = Y[mask]
    else:
        Yaxis = X[mask]
    if len(im.shape) == 3:
        imlen = im[0][mask].shape[0]
    else:
        imlen = im[mask].shape[0]
    keo = np.nan * np.ones((tlen, imlen))
    
    for ii,img in enumerate(im):
        tmp = np.squeeze(img[mask])
        if Xn < 1 and Yn < 1:
            keo[ii] = tmp
        else:
            ix = np.squeeze(np.where(mask == True))[0]
            iy = np.squeeze(np.where(mask == True))[1]
            if Xn >= 1 and Yn < 1:
                direction = 'x'
                N = Xn
            elif Xn < 1 and Yn >= 1:
                direction = 'y'
                N = Yn
            else:
                direction = 'both'
                N = int(np.nanmax(([Xn,Yn])))
            for i, idx in enumerate(ix):
                nbg = getNeighbours(img, idx, iy[i], 
                                    N = N, direction = direction)
                tmp[i] = np.nanmean(nbg)
            keo[ii] = tmp
    
    if filter == 'gaussian':
        keo = ndimage.gaussian_filter(keo, 0.1)
    elif filter =='median':
        keo = ndimage.median_filter(keo, 3)
    else:
        pass
    
    if plotmask:
        # Plot Mask
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax1.pcolormesh(xgrid,ygrid,mask.T)
        ax1.set_xticks(xgrid[::20])
        ax1.grid(color='w', linestyle='--')
        
        ax2 = fig.add_subplot(212)
        ax2.pcolormesh(xgrid,ygrid,mask.T)
        ax2.plot(Y[mask], X[mask], '.r')
        xmin = xgrid[idX].min() - 5
        xmax = xgrid[idX].max() + 5
        ymin = ygrid[idY].min() - 5
        ymax = ygrid[idY].max() + 5
        xstep = np.round((xmax - xmin) / 5)
        ystep = np.round((ymax - ymin) / 5)
        
        ax2.set_xlim([xmin, xmax])
        ax2.set_ylim([ymin, ymax])
        ax2.set_xticks(np.arange(xmin, xmax+1, xstep))
        ax2.set_yticks(np.arange(ymin, ymax+1, ystep))
        ax2.grid(color='w', linestyle='--')
    
    return {'time': time, 'Y': Yaxis, 'keo': keo}

def getTimeSeries(time: Union[list, np.ndarray] = None,
               tlim: Union[list, np.ndarray] = None,
               xgrid: Union[list, np.ndarray] = None,
               ygrid: Union[list, np.ndarray] = None,
               imarray: np.ndarray = None,
               ix: Union[int, float] = None,
               iy: Union[int, float] = None,
               intsize: int = 0):
    if tlim is not None:
        assert type(tlim[0]) == type(time[0])
        idT = (time >= tlim[0]) & (time <= tlim[1])
        time = time[idT]
        imarray = imarray[idT]
    
    if len(xgrid.shape) == 1:
        idX = abs(xgrid - ix).argmin()
    elif len(xgrid.shape) == 2:
        if np.unique(xgrid[0,:]).shape[0] > 1:
            idX = abs(xgrid[0,:] - ix).argmin()
        else:
            idX = abs(xgrid[:,0] - ix).argmin()
    else:
        raise ('Invalid shape for xgrid')
    if len(ygrid.shape) == 1:
        idY = abs(ygrid - iy).argmin()
    elif len(ygrid.shape) == 2:
        if np.unique(ygrid[:,0]).shape[0] > 1:
            idY = abs(ygrid[:,0] - iy).argmin()
        else:
            idY = abs(ygrid[0, 1] - iy).argmin()
    else:
        raise ('Invalid shape for xgrid')
        
    if intsize < 1:
        timeseries = imarray[:,idX,idY]
    else:
        timeseries = np.nan * np.ones(time.shape[0], dtype=np.float16)
        for i in range(time.shape[0]):
            cluster = getNeighbours(imarray[i], idX, idY, N=intsize)
            timeseries[i] = np.nanmean(cluster)
    
    return {'t': time, 'y': timeseries}

def tdft(t,y,T=30,nfft=1024,Nw=240,Nskip=1,window='hamming'):
    Wn = signal.get_window(window,Nw)
    f = np.fft.fftfreq(nfft,d=T) * 1e3 # to mHz
    f = f[1:int(nfft/2)]
    
    Treducted = t[:-Nw]
    Tspecto = Treducted[::Nskip] + timedelta(seconds=Nw/2*T)
    
    for i in np.arange(0,y.shape[0]-Nw,Nskip):
        Stmp = np.fft.fft(y[i:i+Nw]*Wn,n=nfft)
        Sx1 = abs(Stmp[1:int(nfft/2)])**2
        if i == 0:
            Sx = Sx1
        else:
            Sx = np.vstack((Sx,Sx1))
            
    return Tspecto, f, Sx