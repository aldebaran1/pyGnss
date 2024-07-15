#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 16:32:54 2019

@author: smrak
"""

import h5py
import yaml
import os
from numpy import array, logical_and, round, meshgrid, arange, sort, empty, nan_to_num, isfinite

def makeGrid(ylim=[25,50], xlim=[-110,-80], res=0.5):
    """
    Make a grid for an image with a given boundaries and resolution
    """
    xgrid, ygrid = meshgrid(arange(xlim[0],xlim[1],res),arange(ylim[0],ylim[1]+1,res))
    return xgrid, ygrid

def getData(fn):
    fn = h5py.File(fn, 'r')
    rx = fn['data/rx'][:]
    lat = fn['data/table'][:,0]
    lon = fn['data/table'][:,1]
    rxa = [str(l[0].decode('ascii')).lower() for l in rx]
    return lon, lat, array(rxa)

def filterandsave(fn: str = None, ofn: str = None,
                  latlim: list = None,lonlim: list = None, 
                  density: int = None, resolution:int = None):
    """
    """
    if latlim is None:
        latlim = [-90, 90]
    if lonlim is None:
        lonlim = [-180, 180]
    latlim = array(latlim).astype(float)
    lonlim = array(lonlim).astype(float)
    
    if ofn is None:
        ofn = os.path.split(fn)[0] + os.sep
    if os.path.isdir(ofn):
        ofn = os.path.join(os.path.splitext(ofn)[0], os.path.split(fn)[1]) +  '.yaml'
    if os.path.splitext(ofn)[1] != 'yaml':
        ofn = os.path.splitext(ofn)[0] + '.yaml'
    # Get data and locations
    x,y,r = getData(fn)
    # filter - spatial
    idx = (nan_to_num(x) >= lonlim[0]) & (nan_to_num(x) <= lonlim[1]) & (isfinite(x))
    idy = (nan_to_num(y) >= latlim[0]) & (nan_to_num(y) <= latlim[1]) & (isfinite(y))
    ix = logical_and(idx, idy)
    x = x[ix]
    y = y[ix]
    r = r[ix]
    # Relax
    if density is not None:
        xgrid, ygrid = makeGrid(xlim=lonlim, ylim=latlim, res=resolution)
        rx_ix_map = empty((xgrid.shape), dtype='U30')
        idx = array([abs(xgrid[0,:]-x1).argmin() for x1 in x])
        idy = array([abs(ygrid[:,0]-y1).argmin() for y1 in y])
        for i in range(x.size):
            val = rx_ix_map[idy[i], idx[i]]
            if val == '':
                rx_ix_map[idy[i], idx[i]] = str(i)
            else:
                rx_ix_map[idy[i], idx[i]] += ',{}'.format(i)#str(i)
        mask = rx_ix_map == ''
        idrx = []
        for val in rx_ix_map[~mask]:
            if len(val) > 1:
                for i, v in enumerate(val.split(',')):
                    if i == density: break
                    idrx.append(int(v))
            else:
                idrx.append(int(val))
        idrx = sort(idrx)
        x = x[idrx]
        y = y[idrx]
        r = r[idrx]

    outd = [[str(r[i]), str(round(x[i],2)), str(round(y[i],2))] for i in range(r.size)]
    
    # Dump to yaml
    f = open(ofn, 'w')
    yaml.dump({'total': len(outd), 'rx' : outd}, stream=f, default_flow_style=False, allow_unicode=1)
    f.close()
    print ("List successfully saved to ", ofn)
    return

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('rxlist', type = str, help = 'Input data (hdf5)')
    p.add_argument('--ofn', type = str, help = 'Ndestination folder or output file')
    p.add_argument('--lonlim', help = 'Longitude limits', nargs = 2, default=None)
    p.add_argument('--latlim', help = 'Latitude limits', nargs = 2, default=None)
    p.add_argument('--density', type = int, help='Reduce number to rx# per resolution, default=None')
    p.add_argument('-r', '--resolution', type = int, help='Grid resolution to reduce the rx#. Default=1', default=1)
    P = p.parse_args()
    
    filterandsave(fn = P.rxlist, ofn = P.ofn, latlim=P.latlim, lonlim=P.lonlim,
                  resolution=P.resolution, density=P.density)
