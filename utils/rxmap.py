#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 12:09:46 2018

@author: Sebastijan Mrak <smrak@bu.edu>
"""

import matplotlib.pyplot as plt
import h5py
from numpy import array, float16, nanmin, nanmax
import cartopy.crs as ccrs
from cartomap import geogmap as gm
import os
import yaml

def getCoord(fn):
    fn = h5py.File(fn, 'r')
    lat = fn['data/table'][:,0]
    lon = fn['data/table'][:,1]
    
    return lon, lat

def plotMap(fn,lonlim=None,latlim=None, projection='stereo'):
    if os.path.splitext(fn)[1] in ('.h5', '.hdf5'):
        root, fname = os.path.split(fn)
        if latlim is None:
            latlim = [-89.5,89.5]
        if lonlim is None:
            lonlim = [-179.5,180]
#        else:
        latlim = array(latlim, dtype=float)
        lonlim = array(lonlim, dtype=float)
        gm.plotCartoMap(latlim=latlim,lonlim=lonlim, projection='merc')
        lon, lat = getCoord(fn)
    elif os.path.splitext(fn)[1] in ('.yaml', '.yml'):
        root, fname = os.path.split(fn)
        stream = yaml.safe_load(open(fn, 'r'))
        data = array(stream['rx'])
        lon = data[:, 1].astype(float16)
        lat = data[:, 2].astype(float16)
        # Limits
        if lonlim is None:
            
            X = nanmax(lon)
            x = nanmin(lon)
            if X > 150 and x < -150:
                x0 = X - 5 if X-5 < -180 else -180
                x1 = x + 5 if x+5 < 180 else 180
                
            else:
                x0 = x - 5 if x-5 < -180 else -180
                x1 = X + 5 if X+5 < 180 else 180
            lonlim = [round(x0), round(x1)]
        if latlim is None:
            
            y = nanmin(lat)
            Y = nanmax(lat)
            y0 = y-2 if y-2 > -90 else -90
            y1 = Y+2 if Y+2 < 90 else 90
            
            latlim = [round(y0), round(y1)]
        ax = gm.plotCartoMap(latlim=latlim,lonlim=lonlim, projection=projection,
                        states=1, background_color='gray')
    else:
        print ('Wrong file format')
        return
    
    plt.title(fname)
    ax.scatter(lon,lat, marker='.', c='r', s=25, transform=ccrs.PlateCarree())
    print ('Total {} of receivers'.format(lat.size))
    plt.show()
    
if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    
    p.add_argument('hdffile',type=str)
    p.add_argument('-x', '--lonlim', type=float, nargs=2, default = None)
    p.add_argument('-y', '--latlim', type=float, nargs=2, default = None)
    p.add_argument('--proj', type=str, default = 'stereo')
    
    P = p.parse_args()
    plotMap(P.hdffile, P.lonlim, P.latlim, P.proj)