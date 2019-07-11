#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  9 17:24:03 2018

@author: smrak
"""

from six.moves.urllib.parse import urlparse
import ftplib
import subprocess
import os
import numpy as np
from glob import glob
import platform
from datetime import datetime
from dateutil import parser

def unzip(f, timeout=10):
    head, tail = os.path.split(f)
    print ('Unzipping: ', tail)
    if platform.system() == 'Linux':
        try:
            subprocess.call('gzip -d ' + f, shell=True, timeout=timeout)
        except:
            print ('Problems with: ',tail)
    elif platform.system() == 'Windows':
        try:
            subprocess.call('7z x "{}" -o"{}"'.format(f,head), shell=True,timeout=timeout)
            subprocess.call('del "{}"'.format(f), shell=True, timeout=timeout)
        except:
            print ('Problems with: ', tail)
    return

def download(F, rx, filename):
    # Does the file already exists in the destination directory?
    path, tail = os.path.split(filename)
    flist = sorted(glob(path+'/*'))
    fnlist = np.array([os.path.splitext(f)[0] for f in flist])
    if np.isin(os.path.splitext(filename)[0], fnlist):
        print ("{} already exists".format(tail))
    else:
        print ('Downloading to: ', filename)
        try:
            with open(filename, 'wb') as h:
                F.retrbinary('RETR {}'.format(rx), h.write)
        except:
            pass

def getRinexNav(date:str = None,
                odir:str = None,
                const:str = 'gps'):
    """
    year,doy: integer
    db: the name of the database
    odif: final directory to save the rinex files
    """
    # Parse cors url address
    url =  urlparse('ftp://geodesy.noaa.gov/cors/rinex/')
    # Parse date
    try:
        if len(date.split('-')) == 3:
            dt = parser.parse(date)
        elif len(date.split('-')) == 2:
            dt = datetime.strptime(date, '%Y-%j')
        else:
            print ("Wrong date format. Use 'yyyy-mm-dd' or 'yyyy-jjj'")
            exit()
    except Exception as e:
        raise (e)
    year = str(dt.year)
    doy = dt.timetuple().tm_yday
    # Correct spelling to unify the length (char) of the doy in year (DOY)
    if len(str(doy)) == 2:
        doy = '0' + str(doy)
    elif len(str(doy)) == 1:
        doy = '00' + str(doy)
    elif len(str(doy)) == 3:
        doy = str(doy)
    else: 
        print ('Error - Soomething is wrong with your day of year (DOY)')
    
    # Constellation navigation file variable
    sct = {'gps' : 'n',
           'glonass' : 'g',
           'gallileo' : 'e'}
    # Open a connection to the FTP address
    with ftplib.FTP(url[1],'anonymous','guest',timeout=15) as F:
        rpath = url[2] + '/' + year + '/' + doy + '/'
        F.cwd(rpath)
        urlrx = 'brdc' +doy + '0.' + year[-2:]+sct[const]+'.gz'
        print (urlrx)
        try:
            # urlrx must in in a format "nnnDDD0.YYo.xxx"
            download(F, urlrx, odir+urlrx)
        except Exception as e:
            print (e)
            
    unzip(odir+urlrx)
if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('date', type=str, help='Date format YYYY-mm-dd')
    p.add_argument('dir', type=str, help='destination directory')
    p.add_argument('--sct', type=str, help='constellation type',default='gps')
    
    P = p.parse_args()
    # Get file
    getRinexNav(date = P.date, odir = P.dir, const = P.sct)
