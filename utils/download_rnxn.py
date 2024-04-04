#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  9 17:24:03 2018

@author: smrak
"""

import subprocess
import os, gzip, shutil
import numpy as np
from glob import glob
import platform, urllib.request
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
            subprocess.call('7z x "{}" -o"{}"'.format(f, head), shell=True,timeout=timeout)
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
        
def unzip_rm(f, timeout=5):
    # UNZIP
    with gzip.open(f, 'rb') as f_in:
        with open(f[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    if platform.system() == 'Linux':
        subprocess.call('rm -r {}'.format(f), shell=True, timeout=timeout)
    elif platform.system() =='Windows':
        subprocess.call('del "{}"'.format(f), shell=True, timeout=timeout)
    return

def getRinexNav(date:str = None,
                odir:str = None):
    """
    year,doy: integer
    db: the name of the database
    odif: final directory to save the rinex files
    """
    # Parse cors url address
#    url =  urlparse('ftp://geodesy.noaa.gov/cors/rinex/')
    url = 'https://geodesy.noaa.gov/corsdata/rinex/'
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
    Y = year[-2:]
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
    
    ts = (dt - datetime(1980, 1, 6)).total_seconds()
    gpsweek = int(ts / 60 /60 / 24 / 7)
    weekday = (dt.weekday() + 1 ) % 7
    wwwwd = str(gpsweek) + str(weekday)

    # Use HTTPS acces
    urlnav = f"{url}/{year}/{doy}/brdc{doy}0.{Y}n.gz"
    urlsp3 = f"{url}/{year}/{doy}/igs{wwwwd}.sp3.gz"
    navfile = f'{odir}brdc{doy}0.{Y}n.gz'
    sp3file = f'{odir}igs{doy}0.{Y}sp3.gz'
    print (f'Downloading {urlnav}')
    with urllib.request.urlopen(urlnav, timeout=60) as response, open(navfile, 'wb') as out_file:
        data = response.read() # a `bytes` object
        out_file.write(data)
    unzip_rm(navfile)
    print (f'Downloading {urlsp3}')
    with urllib.request.urlopen(urlsp3, timeout=60) as response, open(sp3file, 'wb') as out_file:
        data = response.read() # a `bytes` object
        out_file.write(data)
    unzip_rm(sp3file)
    return
if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('date', type=str, help='Date format YYYY-mm-dd')
    p.add_argument('dir', type=str, help='destination directory')
#    p.add_argument('--type', type=str, help='Navtype type: gps, glonass, gallileo, sp3', default='gps')
    
    P = p.parse_args()
    # Get file
    getRinexNav(date = P.date, odir = P.dir)
