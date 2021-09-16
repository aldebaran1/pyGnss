#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 13:19:38 2017

@author: Sebastijan Mrak <smrak@gmail.com>
"""

from six.moves.urllib.parse import urlparse
import ftplib
import numpy as np
import yaml
from glob import glob
import os
from datetime import datetime
import subprocess
import platform
from dateutil import parser
import urllib
from bs4 import BeautifulSoup

def download_request(urlpath, filename, force=False):
    # Check is destination directory exists?
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
    # Does the file already exists in the destination directory?
    flist = sorted(glob(path+'/*'))
    fnlist = np.array([os.path.splitext(f)[0] for f in flist])
    
    if not np.isin(filename, fnlist):
        # Do you want to override it?
        print ('Downloading file: {}'.format(tail))
        try:
            with urllib.request.urlopen(urlpath, timeout=60) as response, open(filename, 'wb') as out_file:
                data = response.read() # a `bytes` object
                out_file.write(data)
        except:
            pass
        # Else skip the step
    else:
        if force:
            print ('Downloading file: {}'.format(tail))
            try:
                with urllib.request.urlopen(urlpath, timeout=60) as response, open(filename, 'wb') as out_file:
                    data = response.read() # a `bytes` object
                    out_file.write(data)
            except:
                pass
        else:
            print ('{} File already exists'.format(tail))


def download_cddis(F, rx, filename, force=False):
    # Check is destination directory exists?
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
    # Does the file already exists in the destination directory?
    flist = sorted(glob(path+'/*'))
    fnlist = np.array([os.path.splitext(f)[0] for f in flist])
    if not np.isin(filename, fnlist):
        # Do you want to override it?
        print ('Downloading file: {}'.format(tail))
        try:
            F.retrbinary("RETR " + rx, open(filename, 'wb').write)
        except:
            pass
        # Else skip the step
    else:
        if force:
            print ('Downloading file: {}'.format(tail))
            try:
                F.retrbinary("RETR " + rx, open(filename, 'wb').write)
            except:
                pass
        else:
            print ('{} File already exists'.format(tail))

def download(F, rx, filename,force=False):
    def _dl(F,rx):
        try:
            with open(filename, 'wb') as h:
                F.retrbinary('RETR {}'.format(rx), h.write)
        except Exception as e:
            print (e)
        return
    
    # Check is destination directory exists?
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
    # Does the file already exists in the destination directory?
    flist = sorted(glob(path+'/*'))
    fnlist = np.array([os.path.splitext(f)[0] for f in flist])
    if np.isin(filename, fnlist):
        # Do you want to override it?
        if force:
            _dl(F, rx)
        # Else skip the step
        else:
            print ('{} File already exists'.format(tail))
    # Else Download the file
    else:
        print ('Downloading file: {}'.format(tail))
        _dl(F, rx)
        
def getSingleRxUrl(year, doy, F, db, rxn, hr=False):
    d = []
    F.retrlines('LIST', d.append)
    # Find the files
    # CDDIS db
    stations = []
    if db == 'cddis' and not hr:
        if isinstance(rxn, str):
            match = rxn + doy + '0.' + year[-2:] + 'o.Z'
        elif isinstance(rxn, list):
            match = [r + doy + '0.'+year[-2:]+'o.Z' for r in rxn]
            match = np.array(match)
        ds = [line.split()[-1] for line in d]
        ds = np.array(ds)
        idrx = np.where(np.isin(ds,match))[0]
        if idrx.shape[0] > 0:
            stations = ds[idrx]
    elif db == 'cddis' and hr:
        stations = []
        for line in d:
            l = line.split()[-1]
            if rxn in l: 
                stations.append(l)
        stations = np.array(stations)
        
    # CORS db
    elif db == 'cors':
        match = rxn
        ds = [line.split()[-1] for line in d]
        ds = np.array(ds)
        idrx = np.where(np.isin(ds,match))[0]
        if idrx.shape[0] > 0:
            suffix = doy+'0.'+year[-2:]+'d.gz'
            stations = [st+suffix for st in ds[idrx]]
            stations = np.array(stations)
    # EUREF db
    elif db == 'euref':
        if isinstance(rxn, str):
            match = rxn.upper() + doy + '0.'+year[-2:]+'D.Z'
        elif isinstance(rxn, list):
            match = [r.upper() + doy + '0.' + year[-2:] + 'D.Z' for r in rxn]
            match = np.array(match)
        ds = [line.split()[-1] for line in d]
        ds = np.array(ds)
        idrx = np.where(np.isin(ds,match))[0]
        if idrx.shape[0] > 0:
            stations = ds[idrx]
    # UNAVCO db
    elif db == 'unavco' and not hr:
        if isinstance(rxn, str):
            match = rxn + doy + '0.'+year[-2:]+'d.Z'
        elif isinstance(rxn, list):
            match = [r + doy + '0.'+year[-2:]+'d.Z' for r in rxn]
            match = np.array(match)
        ds = [line.split()[-1] for line in d]
        ds = np.array(ds)
        idrx = np.where(np.isin(ds,match))[0]
        if idrx.shape[0] > 0:
            stations = ds[idrx]
    elif db == 'unavco' and hr:
        if isinstance(rxn, str):
            match = rxn + doy + '0.'+year[-2:]+'d.Z'
        elif isinstance(rxn, (list,np.ndarray)):
            match = [r + doy + '0.'+year[-2:]+'d.Z' for r in rxn]
            match = np.array(match)
        ds = [line.split()[-1] for line in d]
        ds = np.array(ds)
        idrx = np.where(np.isin(ds,rxn))[0]
        if idrx.shape[0] > 0:
            suffix = doy+'0.'+year[-2:]+'d.Z'
            stations = [st+suffix for st in ds[idrx]]
            stations = np.array(stations)
    # Return
    return stations

def getStateList(year, doy, F, db, rxn=None, hr=False):
    if isinstance(rxn, str):
        stations = getSingleRxUrl(year,doy,F,db,rxn=rxn, hr=hr)
        print (stations)
    elif isinstance(rxn, (list, np.ndarray)):
        stations = getSingleRxUrl(year,doy,F,db,rxn=rxn, hr=hr)
        print (stations)
    else:
        d = []
        stations = []
        F.retrlines('LIST', d.append)
        if db == 'cddis' and not hr:
            for line in d:
                arg = line.split()[-1]
                if (arg[-2:] == '.Z') or (arg[-2:] == 'ip') or (arg[-2:] == 'gz'):
                    argclober = arg.split('.')
                    if (len(argclober[0]) == 8):
                        stations.append(arg)
        elif db == 'cddis' and hr:
            stations = []
            for line in d:
                arg = line.split()[-1]
                if (arg[-2:] == '.Z') or (arg[-2:] == 'ip') or (arg[-2:] == 'gz'):
                    stations.append(arg)
            stations = np.array(stations)
        elif db == 'cors':
            for line in d:
                arg = line.split()[-1]
                if (len(arg) == 4):
                    try:
                        rx = arg+str(doy)+'0.'+year[-2:]+'d.gz'
                        stations.append(rx)
                    except:
                        pass
        elif db == 'euref':
            for line in d:
                arg = line.split()[-1]
                if (arg[-2:] == '.Z') or (arg[-2:] == 'ip') or (arg[-2:] == 'gz'):
                    argclober = arg.split('.')
                    if (len(argclober[0]) == 8):
                        stations.append(arg)
                        
        elif db == 'unavco':
            for line in d:
                arg = line.split()[-1]
                if not hr:
                    if (arg[-3:] == 'd.Z'):
                        stations.append(arg)
                else:
                    if (len(arg) == 4):
                        rx = arg+str(doy)+'0.'+year[-2:]+'d.Z'
                        stations.append(rx)
    
    return stations

def getRinexObs(date,
                db:str=None,
                odir:str=None,
                rx=None, dllist=None, 
                fix:bool=False,
                hr:bool = False,
                force:bool=False):
    """
    date: dattetime in str format YYYT-mm-dd
    db: the name of the database
    odif: final directory to save the rinex files
    """
    assert db is not None
    assert odir is not None
    
    # Designator
    if platform.system() == 'Linux':
        des = '/'
    elif platform.system() == 'Windows':
        des = '\\'
    # Dictionary with complementary FTP url addresses
    urllist = {#'cddis': 'ftp://cddis.gsfc.nasa.gov/gnss/data/daily/',
               'cddis': 'https://cddis.nasa.gov/archive/gnss/data/daily/',
               'cddishr': 'ftp://ftp.cddis.eosdis.nasa.gov/gps/data/highrate/',
               #'cors':  'ftp://geodesy.noaa.gov/cors/rinex/',
               'cors': 'https://geodesy.noaa.gov/corsdata/rinex/',
               'chain': 'http://chain.physics.unb.ca/data/gps/data/daily/',
               'euref': 'ftp://epncb.oma.be/pub/obs/',
               #'unavco': 'ftp://data-out.unavco.org/pub/rinex/obs/',
               'unavco': 'https://data.unavco.org/archive/gnss/rinex/obs/',
               #'unavcohr': 'ftp://data-out.unavco.org/pub/highrate/1-Hz/rinex/',
               'unavcohr': 'https://data.unavco.org/archive/gnss/highrate/1-Hz/rinex/',
               'ring': 'ftp://bancadati2.gm.ingv.it:2121/OUTGOING/RINEX30/RING/'}
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
    
    # Complete the save directory path
    if not fix:
        if len(str(dt.month)) == 1:
            month = '0' + str(dt.month)
        else:
            month = str(dt.month)
        if len(str(dt.day)) == 1:
            day = '0' + str(dt.day)
        else:
            day = str(dt.day)
        foldername = month + day
        odir += year + des + foldername + des
        if not os.path.exists(odir):
            if not os.path.exists(odir):
                if not os.path.exists(odir):
                    try:
                        subprocess.call('mkdir "{}"'.format(odir), shell=True)
                    except:
                        print ('Cant make the directory')
            try:
                subprocess.call('mkdir "{}"'.format(odir), shell=True)
            except:
                print ('Cant make the directory')
    # Reasign dllist from yaml into rx [list]
    if dllist is not None and isinstance(dllist,str):
        if dllist.endswith('.yaml'):
            stream = yaml.load(open(dllist, 'r'), Loader=yaml.SafeLoader)
            rx = np.array(stream.get('rx'))
            if len(rx.shape) > 1:
                rx = rx[:,0]
        else:
            exit()
    #CDDIS
    if db == 'cddis':
        ftps = ftplib.FTP_TLS(host='gdc.cddis.eosdis.nasa.gov')
        ftps.login(user='anonymous', passwd='sebastijan.mrak@gmail.com')
        ftps.prot_p()
        rpath = 'gnss/data/daily/' + year + '/' + doy + '/' + Y + 'd/'
        ftps.cwd(rpath)
        rxlist = getStateList(year, doy, ftps, db, rxn=rx)
        print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
        for urlrx in rxlist:
            download_cddis(ftps, urlrx, odir+urlrx,force=force)
            
    elif db == 'chain':
        url = f'{urllist[db]}/{year}/{doy}/{Y}d/'
        
        rxlist = []
        with urllib.request.urlopen(url) as response:
            html = response.read().decode('ascii')
            soup = BeautifulSoup(html, 'html.parser')
            for link in soup.find_all('a'):
                if link.get('href') is not None and len(link.get('href')) == 14:
                    rxlist.append(link.get('href')[:4])
                    
        if isinstance(rx, str):
            irx = np.isin(np.asarray(rxlist), rx)
            rxlist = list(np.asarray(rxlist)[irx]) if np.sum(irx) > 0 else None
        
        if rxlist is not None:
            print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
            for rx in rxlist:
                path = f"{url}/{rx}{doy}0.{Y}d.Z"
                print (path)
                ofn = f'{odir}{rx}{doy}0.{Y}d.Z'
                download_request(urlpath=path, filename=ofn, force=force)
                break
        else:
            print ('{} wasnt found'.format(rx))
    
    elif db == 'unavco':
        if hr:
            dbhr = db+'hr'
            url = f'{urllist[dbhr]}/{year}/{doy}/'
        else:
            url = f'{urllist[db]}/{year}/{doy}/'
        rxlist = []
        with urllib.request.urlopen(url) as response:
            html = response.read().decode('ascii')
            soup = BeautifulSoup(html, 'html.parser')
            for link in soup.find_all('a'):
                if link.get('href') is not None and len(link.get('href')[:-1]) == 4:
                    rxlist.append(link.get('href')[:4])
                    
        if isinstance(rx, str):
            irx = np.isin(np.asarray(rxlist), rx)
            rxlist = list(np.asarray(rxlist)[irx]) if np.sum(irx) > 0 else None
        
        if rxlist is not None:
            print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
            for rx in rxlist:
                path = f"{url}/{rx}/{rx}{doy}0.{Y}d.Z"
                print (path)
                ofn = f'{odir}{rx}{doy}0.{Y}d.Z'
                download_request(urlpath=path, filename=ofn, force=force)
                break
        else:
            print ('{} wasnt found'.format(rx))
    
    elif db == 'cors':
        url = f'{urllist[db]}/{year}/{doy}/'
        
        rxlist = []
        with urllib.request.urlopen(url) as response:
            html = response.read().decode('ascii')
            soup = BeautifulSoup(html, 'html.parser')
            for link in soup.find_all('a'):
                if len(link.get('href')[:-1]) == 4:
                    rxlist.append(link.get('href')[:-1])
        if isinstance(rx, str):
            irx = np.isin(np.asarray(rxlist), rx)
            rxlist = list(np.asarray(rxlist)[irx]) if np.sum(irx) > 0 else None
        
        if rxlist is not None:
            print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
            for rx in rxlist:
                path = f"{url}/{rx}/{rx}{doy}0.{Y}d.gz"
                ofn = f'{odir}{rx}{doy}0.{Y}d.gz'
                download_request(urlpath=path, filename=ofn, force=force)
                break
        else:
            print ('{} wasnt found'.format(rx))
    
    else:
        raise('Wrong database')
#    else:
#        # Open a connection to the FTP address
#        with ftplib.FTP(url[1],'anonymous','guest',timeout=45) as F:
#            if db == 'cors':
#                rpath = url[2] + '/' + year + '/' + doy + '/'
#                F.cwd(rpath)
#                # Get the name of all avaliable receivers in the direcotry
#                rxlist = getStateList(year, doy, F, 'cors', rxn=rx)
#                # Download the data
#                print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
#                for urlrx in rxlist:
#                    print (urlrx)
#                    F.cwd(rpath+urlrx[:4]+'/')
#                    download(F, urlrx, odir+urlrx,force=force)
#            elif db == 'euref':
#                rpath = url[2] + '/' + year + '/' + doy + '/'
#                F.cwd(rpath)
#                # Get the name of all avaliable receivers in the direcotry
#                rxlist = getStateList(year, doy, F, db, rxn=rx)
#                # Download the data
#                print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
#                for urlrx in rxlist:
#                    # urlrx must in in a format "nnnDDD0.YYo.xxx"
#                    download(F, urlrx, odir+urlrx,force=force)
#            elif db == 'unavco':
#                rpath = url[2] + '/' + year + '/' + doy + '/'
#                F.cwd(rpath)
#                # Get the name of all avaliable receivers in the direcotry
#                if not hr:
#                    rxlist = getStateList(year, doy, F, db, rxn=rx, hr=hr)
#                    # Download the data
#                    print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
#                    for urlrx in rxlist:
#                        # urlrx must in in a format "nnnDDD0.YYo.xxx"
#                        download(F, urlrx, odir+urlrx,force=force)
#                else:
#                    rxlist = getStateList(year, doy, F, db, rxn=rx, hr=hr)
#                    # Download the data
#                    print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
#                    for urlrx in rxlist:
#                        F.cwd(rpath+urlrx[:4]+'/')
#                        # urlrx must in in a format "nnnDDD0.YYo.xxx"
#                        download(F, urlrx, odir+urlrx, force=force)
            

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('date', help='2017-5-27, or 2017-251', type=str)
    p.add_argument('db',type=str, help='database acronym. Supporting: cddis, \
                   cors, euref')
    p.add_argument('dir',type=str, help='destination directory')
    p.add_argument('-r', '--rx', type=str, help='download a single file for a \
                   given receiver (4 letters)', default=None)
    p.add_argument('-l', '--dllist', type=str, help='A list of receiver names to \
                   be downloaded in yaml cfg format. Look at the example file \
                   dl_list.yaml', default=None)
    p.add_argument('-f', '--force', help='Download anyway, despite the file already exists',
                   action = 'store_true')
    p.add_argument('--highrate', help = 'High rate data if available', 
                   action='store_true')
    p.add_argument('--fixpath', help = 'Fix the odir path?', 
                   action='store_true')
    
    P = p.parse_args()
    if P.db == 'all':
        a = ['cors', 'cddis', 'euref', 'unavco']
        for db in a:
            getRinexObs(date = P.date, db = db, 
                        odir = P.dir, rx = P.rx, dllist = P.dllist, 
                        hr = P.highrate, force = P.force, fix = P.fixpath)
    elif P.db == 'conus':
        a = ['cors', 'unavco', 'cddis']
        for db in a:
            getRinexObs(date = P.date, db = db, 
                        odir = P.dir, rx = P.rx, dllist = P.dllist, 
                        hr = P.highrate, force = P.force, fix = P.fixpath)
    elif P.db == 'south':
        a = ['cddis', 'unavco']
        for db in a:
            getRinexObs(date = P.date, db = db, 
                        odir = P.dir, rx = P.rx, dllist = P.dllist, 
                        hr = P.highrate, force = P.force, fix = P.fixpath)
    else:
        getRinexObs(date = P.date, db = P.db, odir = P.dir, 
                    rx = P.rx, dllist = P.dllist, hr = P.highrate, 
                    force = P.force, fix = P.fixpath)
            
