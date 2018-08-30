# -*- coding: utf-8 -*-
"""
Created on Sat Aug 25 15:09:07 2018

@author: smrak
"""

import georinex as gr
from datetime import datetime
from pandas import Timestamp
from numpy import array
import xarray

mahfolder = 'E:\\mahali\\2015-10-07\\'
obs = 'mah82800.15o'
nc = '\\nc\\mah82800.nc'
fno = mahfolder+obs
fn = mahfolder+nc

# Open with Xarray
data = xarray.open_dataset(fn, group='OBS')
L1 = data.sel(sv='G23')['L1'].values[1000:]
t64 = data.time.values
dt = array([Timestamp(t).to_pydatetime() for t in t64])

rx_xyz = xarray.open_dataset(fn, group='POSITION').XYZ.values

# Convert to nc
#save = 'E:\\mahali\\2015-10-07\\nc\\'
#data = gr.load(fno, ofn=save)

#
#F = netCDF4.Dataset(fn, 'a')
#fgr = F.createGroup('POSITION')
#fgr.createDimension('ecef', 3)
#xyz = fgr.createVariable('XYZ', 'f4', 'ecef')
#xyz[:] = [100.25, 200.5, 300.0]
#F.close()
#
## Read
#F = netCDF4.Dataset(fn, 'r')
#pos = F.groups['POSITION']
#print (pos.variables['XYZ'][:])
#F.close()
