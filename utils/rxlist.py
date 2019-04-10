#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 12:08:38 2018

@author: Sebastijan Mrak <smrak@bu.edu>
"""
import h5py
import yaml
import os
import glob
import georinex as grx
from numpy import nan, zeros
from pymap3d import ecef2geodetic

def getRxList(folder, sufix):
    sx = sufix.split(',')
    flist = []
    for sx in sx:
        filestr = os.path.join(folder, sx)
        flist += sorted(glob.glob(filestr))
    return flist

def writeRxlist2HDF(obsfolder='/media/smrak/Eclipse2017/Eclipse/cors/all/233/',
                    sufix='*.*d',listfilename=None):
    """
    Make a list of receivers in a given folder, ordered by their geographical
    location. The list is organized as hdf5 file.
    "obsfolder" is path to the folder with RINEX files, which need to be first
    converted into .yaml cfg file containing their header information. 
    "listfilename" is the wanted filename which will contain converted lat/lon
    information. Itshould contain path/to/file.h5, if you forget .h5 extension,
    it will auto make one for you.
    """
    rxlist = getRxList(obsfolder, sufix)
    if listfilename is None:
        listfilename = obsfolder
    assert len(rxlist) > 0
    head, tail = os.path.split(listfilename)
    year = tail[-3:-1]
    doy = tail[4:7]
    if tail == '':
        F = os.path.expanduser(listfilename).split(os.sep)
        year = F[-3]
        doy = F[-2]
        listfilename = listfilename + 'rxlist{}.{}.h5'.format(doy,year)
    if listfilename[-3:] != '.h5':
        listfilename += '.h5'
    print ('Number of receivers in the folder: ', len(rxlist))
    c = 0
    table = nan*zeros((len(rxlist),2))
    for fn in rxlist:
        sx = os.path.splitext(fn)[1]
        try:
            if sx.endswith('.yaml'):
                stream = yaml.load(open(fn, 'r'))
                rx_xyz = stream.get('APPROX POSITION XYZ')
                rec_lat, rec_lon, rec_alt = ecef2geodetic(rx_xyz[0], rx_xyz[1], rx_xyz[2])
            
            elif sx.endswith('o') or sx.endswith('d'):
                hdr = grx.rinexheader(fn)
                if 'position_geodetic' in hdr:
                    rec_lat, rec_lon, rec_alt = hdr['position_geodetic']
                elif 'position' in hdr:
                    rx_xyz = hdr['position']
                    rec_lat, rec_lon, rec_alt = ecef2geodetic(rx_xyz[0], rx_xyz[1], rx_xyz[2])
                elif 'APPROX POSITION XYZ':
                    rx_xyz = hdr['APPROX POSITION XYZ']
                    rec_lat, rec_lon, rec_alt = ecef2geodetic(rx_xyz[0], rx_xyz[1], rx_xyz[2])
                else:
                    print ('Cant find the position of the gps receiver: {}'.format(fn))
                    rec_lat = nan
                    rec_lon = nan
            table[c,:] = [rec_lat, rec_lon]
            c+=1
        except Exception as e:
            print ('Couldnt process: ', fn)
            c+=1
    print ('Saving as: ', listfilename)
    gpslist = [os.path.split(l)[1][:4] for l in rxlist]
    h5file = h5py.File(listfilename, 'w')
    tab = h5file.create_group('data')
    tab.create_dataset('table', data=table)
    asciiList = [n.encode("ascii", "ignore") for n in gpslist]
    tab.create_dataset('rx', (len(asciiList),1),'S10', asciiList)
    h5file.close()
    print ('Successfully saved!')

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('folder',type=str)
    p.add_argument('-o', '--odir',type=str, help='Output directory', default=None)
    p.add_argument('-n', '--filename', help='list filename',type=str, default='')
    p.add_argument('-s', '--sufix', 
                   help='obs suffix to take? *.*o; *.*d; or *.yaml (default = *.*d,*.*o)',
                   type=str, default='*.*d,*.*o')
    P = p.parse_args()
    

    if P.filename == '':
        writeRxlist2HDF(obsfolder=P.folder, listfilename=P.odir, sufix=P.sufix)
    else: 
        writeRxlist2HDF(obsfolder=P.folder, listfilename=P.odir+P.filename, sufix=P.sufix)