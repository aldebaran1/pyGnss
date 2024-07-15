#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 13:25:18 2024

@author: mraks1
"""
import glob
import xarray as xr
import georinex as gr
import numpy as np 
from datetime import datetime

date = datetime(2023, 12, 16)
yyyy = date.year
yy = date.strftime("%y")
doy = date.strftime("%j")
rx = 'pove'

root = '/Users/mraks1/Library/CloudStorage/OneDrive-UCB-O365/Projects/AFOSR-AERO-2024/rinex/'
rxo_hh = sorted(glob.glob(root + f'{rx.upper()}00BRA_R_{yyyy}{doy}0000_01D_01S_MO.crx/gnss/data/highrate/{yyyy}/{doy}/{yy}d/*'))
for ihh, hh in enumerate(rxo_hh):
    if ihh == 0:
        rxoo = np.array(sorted(glob.glob(hh+'/*.gz')))
    else:
        rxoo = np.hstack((rxoo, np.array(sorted(glob.glob(hh+'/*.gz')))))
        

for ioo, f in enumerate(rxoo):
    print (f"{ioo+1}/{rxoo.size}")
    if ioo == 0:
        X = gr.load(rxoo[ioo], use='G', fast=True, )
    else:
        X = xr.concat((X, gr.load(rxoo[ioo], use='G', fast=True)), dim='time', data_vars=('L1C', 'C1C', 'L2W', 'C2W'))

X['filename'] = X.filename.replace("15M", "01D")
ofn = root + X.filename.values.item() + '.nc'
X.to_netcdf(ofn)