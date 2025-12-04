#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 13:28:09 2019

@author: smrak
"""

from six.moves.urllib.parse import urlparse
import ftplib
from glob import glob
import numpy as np
import subprocess
import os
import platform
from datetime import datetime, timedelta
from dateutil import parser

def unzip(f, timeout=10):
    head, tail = os.path.split(f)
    print ('Unzipping: ', tail)
    if platform.system() in ('Linux', 'Darwin'):
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
    
    if not os.path.exists(path):
        #NO? Well, create the directory. 
        try:
            if platform.system() == 'Linux':
                subprocess.call('mkdir -p {}'.format(path), shell=True)
            elif platform.system() == 'Windows':
                subprocess.call('mkdir "{}"'.format(path), shell=True)
        except:
            print ('Cant make the directory')

    print ('Downloading file: {}'.format(tail))
    try:
        F.retrbinary("RETR " + rx, open(filename, 'wb').write)
    except:
        pass


def getRinexNav(date:str = None,
                odir:str = None):
    """
    year,doy: integer
    db: the name of the database
    odif: final directory to save the rinex files
    """
    # Parse cors url address
    # url =  urlparse('ftp://cddis.nasa.gov/gnss/products/ionex/')
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
    yy = dt.strftime("%y")
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
        
        
    ftps = ftplib.FTP_TLS(host='gdc.cddis.eosdis.nasa.gov')
    # ftps.connect(port=43094)
    ftps.login(user='anonymous', passwd='sebastijan.mrak@gmail.com')
    ftps.prot_p()
    rpath = 'gnss/products/ionex/' + year + '/' + doy + '/'
    ftps.cwd(rpath)
    urlrx_new = 'ESA0OPSRAP_' + year + doy + '0000_01D_01H_GIM.INX.gz' 
    urlrx_old = 'esag' + doy + '0.' + yy + 'i.Z'
    d = []
    ftps.retrlines('LIST', d.append)
    files = [a.split()[-1] for a in d] 
    if urlrx_new in files:
        download(ftps, urlrx_new, odir+urlrx_new)
        unzip(odir+urlrx_new)
    elif urlrx_old in files:
        download(ftps, urlrx_old, odir+urlrx_old)
        unzip(odir+urlrx_old)
    else:
        pass
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
