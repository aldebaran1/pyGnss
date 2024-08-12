#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  9 17:24:03 2018

@author: smrak
"""

import subprocess
import ftplib
import os, gzip, shutil
import numpy as np
from glob import glob
import platform, urllib.request
from datetime import datetime, timedelta
from dateutil import parser
import ssl

ssl._create_default_https_context = ssl._create_unverified_context

def unzip(f, timeout=10):
    head, tail = os.path.split(f)
    print ('Unzipping: ', tail)
    if platform.system() == 'Linux':
        try:
            subprocess.call('gzip -d -f -q' + f, shell=True, timeout=timeout)
        except:
            print ('Problems with: ',tail)
    elif platform.system() == 'Windows':
        try:
            subprocess.call('7z x "{}" -o"{}"'.format(f, head), shell=True, timeout=timeout)
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
    wwww = int(ts / 60 /60 / 24 / 7)
    weekday = (dt.weekday() + 1 ) % 7
    wwwwd = str(wwww) + str(weekday)
    
    # Rinex3 Mixed Navigation File"
    ftps = ftplib.FTP_TLS(host='gdc.cddis.eosdis.nasa.gov')
    ftps.login(user='anonymous', passwd='sebastijan.mrak@gmail.com')
    ftps.prot_p()
    parent_dir = ftps.pwd()
    rpath = f'gnss/data/daily/{year}/brdc/'
    ftps.cwd(rpath)
    target = np.array([f'brdc{doy}0.{Y}n.Z', f'BRDC00IGS_R_{year}{doy}0000_01D_MN.rnx.gz'])
    d = []
    ftps.retrlines('LIST', d.append)
    dta = np.array([dd.split()[-1] for dd in d])
    isin = np.isin(dta, target)
    if np.sum(isin) > 0:
        for f in dta[isin]:
            if not os.path.exists(odir):
                subprocess.call(f'mkdir -p "{odir}"', shell=True)
            print (f"Downloading: {rpath}{f}")
            ftps.retrbinary("RETR " + f, open(f'{odir}{f}', 'wb').write)
            unzip_rm(f'{odir}{f}')
    
    # Rinex3 Mixed Orbits in SP3 format
    ftps.cwd(parent_dir)
    rpath = f'gnss/products/{wwww}/'
    ftps.cwd(rpath)
    target = f'GFZ0MGXRAP_{year}{doy}0000_01D_05M_ORB.SP3.gz'
    d = []
    ftps.retrlines('LIST', d.append)
    dta = np.array([dd.split()[-1] for dd in d])
    isin = np.isin(dta, target)
    if np.sum(isin) > 0:
        for f in dta[isin]:
            if not os.path.exists(odir):
                subprocess.call(f'mkdir -p "{odir}"', shell=True)
            print (f"Downloading: {rpath}{f}")
            ftps.retrbinary("RETR " + f, open(f'{odir}{f}', 'wb').write)
            unzip_rm(f'{odir}{f}')
    

    # Use HTTPS acces
    urlnav = f"{url}/{year}/{doy}/brdc{doy}0.{Y}n.gz"
    urlsp3 = f"{url}/{year}/{doy}/igs{wwwwd}.sp3.gz"
    # urlsp3_r3 = f'{url}/{year}/{doy}/IGS0OPSRAP_{dt.strftime("%Y%j")}0000_01D_15M_ORB.SP3.gz'
    navfile = f'{odir}brdc{doy}0.{Y}n.gz'
    sp3file = f'{odir}igs{doy}0.{Y}sp3.gz'
    # sp3file_r3 = f'{odir}IGS0OPSRAP_{dt.strftime("%Y%j")}0000_01D_15M_ORB.SP3.gz'
    
    print (f'Downloading {urlnav}')
    with urllib.request.urlopen(urlnav, timeout=60) as response, open(navfile, 'wb') as out_file:
        if not os.path.exists(odir):
            subprocess.call(f'mkdir -p "{odir}"', shell=True)
        data = response.read() # a `bytes` object
        out_file.write(data)
    unzip_rm(navfile)
    print (f'Downloading {urlsp3}')
    with urllib.request.urlopen(urlsp3, timeout=60) as response, open(sp3file, 'wb') as out_file:
        if not os.path.exists(odir):
            subprocess.call(f'mkdir -p "{odir}"', shell=True)
        data = response.read() # a `bytes` object
        out_file.write(data)
    unzip_rm(sp3file)
    
    return

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('date', type=str, help='Date format YYYY-mm-dd in comma-separated start,end if wanted')
    p.add_argument('dir', type=str, help='destination directory')
    
    P = p.parse_args()
    dates = P.date.split(',')
    if len(dates) < 1:
        print ("the date is a required argument")
        exit()
    elif len(dates) == 2:
        dtdates = np.arange(parser.parse(dates[0]),parser.parse(dates[1])+timedelta(hours=1), timedelta(days=1)).astype('datetime64[s]').astype(datetime)
        dates = [t.strftime("%Y-%m-%d") for t in dtdates]
    elif len(dates) > 2:
        print ("Can't compile the dates")
        exit()
    for d in dates:
        # Get file
        getRinexNav(date = d, odir = P.dir)
