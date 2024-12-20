#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 13:19:38 2017

@author: Sebastijan Mrak <smrak@gmail.com>
"""

# from six.moves.urllib.parse import urlparse
import requests
import ftplib
import numpy as np
import yaml
from glob import glob
import os
from datetime import datetime, timedelta
import subprocess
import platform
from dateutil import parser
import urllib.request
from bs4 import BeautifulSoup
from earthscope_sdk.auth.device_code_flow import DeviceCodeFlowSimple
from earthscope_sdk.auth.auth_flow import NoTokensError
from pathlib import Path
import ssl
import string
import warnings

ssl._create_default_https_context = ssl._create_unverified_context
#Change to your preference

warnings.filterwarnings("ignore")

token_path = os.path.expanduser("~") + '/pyGnss/utils/'
token_path= os.getcwd() + os.sep

hhindd = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23']

urllist = {'cddis': 'https://cddis.nasa.gov/archive/gnss/data/daily/',
           'cddishr': 'https://cddis.nasa.gov/archive/gnss/data/highrate/',
           'cors': 'https://geodesy.noaa.gov/corsdata/rinex/',
           'chain': 'http://chain.physics.unb.ca/data/gps/data/daily/',
           'chainhr': 'http://chain.physics.unb.ca/data/gps/data/highrate/',
           'euref': 'https://epncb.oma.be/pub/RINEX/',
           'eurefhr': 'https://igs.bkg.bund.de/root_ftp/EUREF/highrate/',
           'chile': 'https://gps.csn.uchile.cl/',
           'brasil': 'https://geoftp.ibge.gov.br/informacoes_sobre_posicionamento_geodesico/rbmc/dados/',
           'unavco': 'https://data.unavco.org/archive/gnss/rinex/obs/',
           'unavcohr': 'https://data.unavco.org/archive/gnss/highrate/1-Hz/rinex/',
           'ring': 'bancadati2.gm.ingv.it', #italy
           'bev': 'https://gnss.bev.gv.at/at.gv.bev.dc/data/', #austria
           'sonel': 'ftp.sonel.org',
           'nz': 'https://data.geonet.org.nz/gnss/rinex/',
           'nl': 'https://gnss1.tudelft.nl/dpga/rinex/',
           'd': 'https://igs.bkg.bund.de/root_ftp/GREF/obs/',
           'fr': 'rgpdata.ensg.eu/data/',
           'fr2': 'ftp://renag.unice.fr/',
           'es': 'https://datos-geodesia.ign.es/ERGNSS/diario_30s/',
           'epos': 'https://datacenter.gnss-epos.eu/',
           's': 'ftpswepos-open.lantmateriet.se'
           }
    
def download_request(urlpath, filename, force=False, hr=False):
    # Check is destination directory exists?
    path, tail = os.path.split(filename)

    if not os.path.exists(path):
        #NO? Well, create the directory. 
        try:
            if platform.system() in ('Linux', 'Darwin'):
                subprocess.call('mkdir -p "{}"'.format(path), shell=True)
            elif platform.system() == 'Windows':
                subprocess.call('mkdir "{}"'.format(path), shell=True)
            else:
                print ("Unknown system; terminating the operation")
                exit()
        except:
            print ('Cant create the directory')
    # Does the file already exists in the destination directory?
    flist = sorted(glob(path+'/*'))
    fnlist = np.array([os.path.splitext(f)[0] for f in flist])
    if not hr:
        fnamelist = np.array([os.path.split(f)[-1][:4].lower() for f in fnlist])
        fn = os.path.split(filename)[-1][:4].lower()
    else:
        fn = filename
        fnamelist = fnlist
    if not np.isin(fn, fnamelist):
        # Do you want to override it?
        print ('Downloading file: {}'.format(tail))
        try:
            with urllib.request.urlopen(urlpath, timeout=60) as response, open(filename, 'wb') as out_file:
                data = response.read() # a `bytes` object
                out_file.write(data)
        except Exception as e:
            print (e)
        # Else skip the step
    else:
        if force:
            print ('Downloading file: {}'.format(tail))
            try:
                with urllib.request.urlopen(urlpath, timeout=60) as response, open(filename, 'wb') as out_file:
                    data = response.read() # a `bytes` object
                    out_file.write(data)
            except Exception as e:
                print (e)
        else:
            print ('{} File already exists'.format(tail))


def download_cddis(F, rx, filename, force=False):
    # Check is destination directory exists?
    path, tail = os.path.split(filename)
    if not os.path.exists(path):
        #NO? Well, create the directory. 
        try:
            if platform.system() in ('Linux', 'Darwin'):
                subprocess.call('mkdir -p {}'.format(path), shell=True)
            elif platform.system() == 'Windows':
                subprocess.call('mkdir "{}"'.format(path), shell=True)
        except:
            print ('Cant create the directory')
    # Does the file already exists in the destination directory?
    flist = sorted(glob(path+f'{os.sep}*'))
    fnlist = np.array([os.path.split(f)[1][:4].lower() for f in flist])
    if not np.isin(os.path.split(filename)[1][:4].lower(), fnlist):
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
                    if (len(argclober[0]) == 8) or (len(argclober[0].split('_'))==6):
                        stations.append(arg)
                        
        elif db == 'cddis' and hr:
            stations = []
            if len(d) < 20 and len(d[0].split()[-1]) == 3:
                F.cwd(f'{year[-2:]}d/')
                for i, hh in enumerate(hhindd):
                    # print (i, hh)
                    d = []
                    F.cwd(f'{hh}/')
                    F.retrlines('LIST', d.append)
                    # print (d)
                    for line in d:
                        arg = line.split()[-1]
                        if arg[-3:] == '.gz':
                            stations.append(f"{year[-2:]}d/{hh}/{arg}")
                            # stations.append(f'{hh}/{arg}')
                    F.cwd('../')
                F.cwd('../')
                    
            
            else:
                for line in d:
                    arg = line.split()[-1]
                    if arg[-3:] == 'tar':
                        if "MO.crx.tar" in arg.split('_'):
                            stations.append(arg)
                        elif f"{year[-2:]}d" in arg.split('.'):
                            stations.append(arg)
                        else:
                            pass
        elif db == 'cors':
            for line in d:
                arg = line.split()[-1]
                if (len(arg) == 4):
                    try:
                        rx = arg+str(doy)+'0.'+year[-2:]+'d.gz'
                        stations.append(rx)
                    except:
                        pass
        elif db == 'ring':
            for line in d:
                arg = line.split()[-1]
                if arg[-2:] in('gz', '.Z', 'ip'):
                    stations.append(arg)
            stations = np.array(stations)
        elif db == 'euref':
            for line in d:
                arg = line.split()[-1]
                if (arg[-2:] == '.Z') or (arg[-2:] == 'ip') or (arg[-2:] == 'gz'):
                    argclober = arg.split('.')
                    if (len(argclober[0]) == 8):
                        stations.append(arg)
        elif db == 'sonel':
            for line in d:
                arg = line.split()[-1]
                if arg.endswith(('d.Z', 'MO.crx.gz')):
                    stations.append(arg)
        
        elif db in ('fr', 's'):
            for line in d:
                arg = line.split()[-1]
                if arg.endswith( '01D_30S_MO.crx.gz'):
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
    des = os.sep
    # Dictionary with complementary FTP url addresses

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
            if platform.system() == "Windows":
                subprocess.call('mkdir "{}"'.format(odir), shell=True)
            else:
                subprocess.call('mkdir -p "{}"'.format(odir), shell=True)
            # try:
                # subprocess.call('mkdir "{}"'.format(odir), shell=True)
            # except:
                # print ('Cant make the directory')
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
        #This one is tricky because it contains Rinex2 and Rinex3 files that do not overlap
        if not hr:
            rpath = 'gnss/data/daily/' + year + '/' + doy + '/' + Y + 'd/'
        else:
            rpath = 'gnss/data/highrate/' + year + '/' + doy + '/'
        ftps.cwd(rpath)
        rxlist = np.array(getStateList(year, doy, ftps, db, hr=hr))

        
        if not hr:
            rxnames = np.array([r[:4].lower() for r in rxlist])
            _, ix = np.unique(rxnames, return_index=True)
            rxlist = rxlist[ix]
            rxnames = rxnames[ix]
        else:
            rxnames = np.array([r.split("/")[-1][:4].lower() for r in rxlist])
            
        # print (rxnames)
        if isinstance(rx, str):
            irx = np.isin(rxnames, rx.lower())
            rxlist = rxlist[irx]
            rxlen = np.array([len(r) for r in rxlist])
            lens, _ = np.unique(rxlen, return_index=True)
            if lens.size > 1:
                keeper = max(lens)
                ilens = np.isin(rxlen, keeper)
                rxlist = rxlist[ilens]
        print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
        if not hr:
            for urlrx in rxlist:
                download_cddis(ftps, urlrx, odir+urlrx,force=force)
        else:
            for rx in rxlist:
                pth = rx.split('/')
                ftps.cwd(f'{pth[0]}/{pth[1]}/')
                download_cddis(ftps, pth[2], odir+pth[2], force=force)
                ftps.cwd('../')
                ftps.cwd('../')
            
    elif db == 'chain':
        rxlist = []
        if hr:
            rxpath = []
            db += 'hr'
            url = f'{urllist[db]}/{year}/{doy}/'
            for i,hh in enumerate(hhindd):
                urla = url + hh + '/'
                with urllib.request.urlopen(urla) as response:
                    html = response.read().decode('ascii')
                    soup = BeautifulSoup(html, 'html.parser')
                    for link in soup.find_all('a'):
                        if link.get('href') is not None and len(link.get('href')) == 16:
                            rxpath.append(f'{hh}/')
                            rxlist.append(link.get('href'))
            rxlist = np.array(rxlist)
            rxpath = np.array(rxpath)
            rxnames = np.array([r[:4] for r in rxlist])
            
            if isinstance(rx, str):
                irx = np.isin(rxnames, rx)
                rxlist = rxlist[irx] if np.sum(irx) > 0 else None
                rxpath = rxpath[irx] if np.sum(irx) > 0 else None
        else:
            url = f'{urllist[db]}/{year}/{doy}/{Y}d/'
            with urllib.request.urlopen(url) as response:
                html = response.read().decode('ascii')
                soup = BeautifulSoup(html, 'html.parser')
                for link in soup.find_all('a'):
                    if link.get('href') is not None and len(link.get('href')) == 14:
                        rxlist.append(link.get('href'))
            rxlist = np.array(rxlist)
            rxnames = np.array([r[:4] for r in rxlist])
            if isinstance(rx, str):
                irx = np.isin(rxnames, rx)
                rxlist = rxlist[irx] if np.sum(irx) > 0 else None
        
        if rxlist is not None:
            print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
            for i, rx in enumerate(rxlist):
                if hr:
                    path = f"{url}/{rxpath[i]}/{rx}"
                else:
                    path = f"{url}/{rx}"
                ofn = f"{odir}{rx}"    
                download_request(urlpath=path, filename=ofn, force=force, hr=hr)
        else:
            print ('{} wasnt found'.format(rx))
    
    elif db == 'unavco':
        device_flow = DeviceCodeFlowSimple(Path(token_path))
        try:
            # get access token from local path
            device_flow.get_access_token_refresh_if_necessary()
        except NoTokensError:
            # if no token was found locally, do the device code flow
            device_flow.do_flow()
        token = device_flow.access_token
        
        rxlist = []
        if hr:
            db = db+'hr'
            url = f'{urllist[db]}{year}/{doy}/'
            r = requests.get(url, headers={"authorization": f"Bearer {token}"},verify=False)
            d = []
            if r.status_code == requests.codes.ok:
                for data in r:
                    d.append(data.decode('ascii').replace('\n',''))
            soup = BeautifulSoup("".join(d), 'html.parser')
            for link in soup.find_all('a'):
                if len(link.get('href').split('/')[-2]) == 4:
                    rxlist.append(link.get('href').split('/')[-2])
            rxlist = np.unique(rxlist)
            if isinstance(rx, str):
                irx = np.isin(np.asarray(rxlist), rx)
                rxlist = list(np.asarray(rxlist)[irx]) if np.sum(irx) > 0 else None
            
            if rxlist is not None:
                print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
                flist = sorted(glob(odir+os.sep+'*'))
                fnames = np.array([os.path.split(f)[1][:4] for f in flist])
                for rx in rxlist:
                    if np.isin(rx, fnames) and (not force):
                        print ('{} already exists'.format(rx))
                        continue
                    path = f"{url}/{rx}/{rx}{doy}0.{Y}d.Z"
                    ofn = f'{odir}{rx}{doy}0.{Y}d.Z'
                    if not os.path.exists(odir):
                        print ("Making new directory:\n", odir)
                        if platform.system() in ('Linux', 'Darwin'):
                            subprocess.call(f'mkdir -p "{odir}"', shell=True)
                        else:
                            subprocess.call(f'mkdir "{odir}"', shell=True)
                    with open(ofn, 'wb') as f:
                        print (f"Downloading {path}:")
                        r = requests.get(path, headers={"authorization": f"Bearer {token}"}, verify=False)
                        for data in r:
                            f.write(data)
                    f.close()
            else:
                print ('{} wasnt found'.format(rx))
        else:
            url = f'{urllist[db]}/{year}/{doy}/'
        
            r = requests.get(url, headers={"authorization": f"Bearer {token}"}, verify=False)
            d = []
            if r.status_code == requests.codes.ok:
                for data in r:
                    d.append(data.decode('ascii').replace('\n',''))
            soup = BeautifulSoup("".join(d), 'html.parser')
            for link in soup.find_all('a'):
                if link.get('href') is not None and len(link.get('href').split('/')[-1]) == 14:
                    rxlist.append(link.get('href').split('/')[-1][:4])
            rxlist = np.unique(rxlist)
            if isinstance(rx, str):
                irx = np.isin(np.asarray(rxlist), rx)
                rxlist = list(np.asarray(rxlist)[irx]) if np.sum(irx) > 0 else None
            
            if rxlist is not None:
                print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
                flist = sorted(glob(odir+os.sep+'*'))
                fnames = np.array([os.path.split(f)[1][:4] for f in flist])
                for rx in rxlist:
                    if np.isin(rx, fnames) and (not force):
                        print ('{} already exists'.format(rx))
                        continue
                    if hr:
                        path = f"{url}/{rx}/{rx}{doy}0.{Y}d.Z"
                        ofn = f'{odir}{rx}{doy}0.{Y}d.Z'
                    else:
                        path = f"{url}/{rx}{doy}0.{Y}d.Z"
                        ofn = f'{odir}{rx}{doy}0.{Y}d.Z'
                    if not os.path.exists(odir):
                        print ("Making new directory:\n", odir)
                        if platform.system() in ('Linux', 'Darwin'):
                            subprocess.call(f'mkdir -p "{odir}"', shell=True)
                        else:
                            subprocess.call(f'mkdir "{odir}"', shell=True)
                    with open(ofn, 'wb') as f:
                        print (f"Downloading {path}:")
                        r = requests.get(path, headers={"authorization": f"Bearer {token}"}, verify=False)
                        for data in r:
                            f.write(data)
                    f.close()
            else:
                print ('{} wasnt found'.format(rx))

    elif db == 'cors':
        url = f'{urllist[db]}/{year}/{doy}/'
        if hr:
            print ("CORS does not support highrate data files")
            return 
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
                download_request(urlpath=path, filename=ofn, force=force, hr=hr)
        else:
            print ('{} wasnt found'.format(rx))
    
    elif db == 'chile':
        """
        This thing doesn't work anymore becuase the wbpage uses javascript to generate the webpage
        """
        # url = f'{urllist[db]}data/{year}/{doy}/'
        # stations_url = f'{urllist[db]}/CSN_GNSS.info'
        # print (url)
        # headers = {
        #     'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.36',
        #     'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
        #     'Accept-Language': 'en-US,en;q=0.5',
        #     'Connection': 'keep-alive',
        #     'Upgrade-Insecure-Requests': '1',
        # }
        # req = urllib.request.Request(stations_url, headers=header)
        # s = requests.Session()
        # s.get(urllist[db], timeout=60, headers=headers, verify=False)
        # s.raise_for_status()
        # print (s.get(urllist[db], timeout=60, headers=headers, verify=False))
        # s.get(url, timeout=60, headers=headers, verify=False)
        # print (s.get(url, timeout=60, headers=headers, verify=False))
        return
        # with urllib.request.urlopen(req, timeout=60) as response:
        #     html = response.read().decode('ascii').split('\n')
        #     rxlist = []
        #     for i,l in enumerate(html):
        #         if i < 2:
        #             continue
        #         try:
        #             rxn = l[:4]
        #             if len(rxn) == 4:
        #                 rxlist.append(rxn.lower())
        #         except:
        #             pass
        #         del rxn
        # rxlist = np.unique(np.asarray(rxlist))
        # if isinstance(rx, str):
        #     irx = np.isin(rxlist, rx)
        #     rxlist = rxlist[irx] if np.sum(irx) > 0 else None
        
        # if rxlist is not None:
        #     print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
        #     for rx in rxlist:
        #         path = f"{url}/{rx}{doy}0.{Y}d.Z"
        #         ofn = f'{odir}{rx}{doy}0.{Y}d.Z'
        #         download_request(urlpath=path, filename=ofn, force=force, hr=hr)
        # else:
        #     print ('{} wasnt found'.format(rx))
            
    elif db == 'euref':
        if hr:
            db+='hr'
        
        url = f'{urllist[db]}/{year}/{doy}/'
        
        if hr:
            rxlist = []
            aainhh = np.array(list(string.ascii_lowercase[:23]))
            for i, aa in enumerate(aainhh): 
                urla = url + f'{aa}/'
                with urllib.request.urlopen(urla) as response:
                    html = response.read().decode('ascii')
                    soup = BeautifulSoup(html, 'html.parser')
                    for link in soup.find_all('a'):
                        if len(link.get('href').split('.')[0].split("_")) == 6:
                            rxlist.append(f"{aa}/{link.get('href')}")
            rxlist = np.array(rxlist)
            rxnames = np.array([r.split('/')[-1][:4].lower() for r in rxlist])
            
            if isinstance(rx, str):
                irx = np.isin(rxnames, rx.lower())
                rxlist = rxlist[irx] if np.sum(irx) > 0 else None
                        
        else:
            with urllib.request.urlopen(url) as response:
                html = response.read().decode('ascii')
                soup = BeautifulSoup(html, 'html.parser')
                rxlist = []
                for link in soup.find_all('a'):
                    if len(link.get('href').split('.')[0].split("_")) == 6:
                        rxlist.append(link.get('href'))
            rxlist = np.array(rxlist)
            rxnames = np.array([r[:4].lower() for r in rxlist])
                                
            if isinstance(rx, str):
                irx = np.isin(rxnames, rx.lower())
                rxlist = rxlist[irx] if np.sum(irx) > 0 else None
                
        #    Download the data
        print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
        for rx in rxlist:
            if not hr:
                ofn = f'{odir}/{rx}'
            else:
                ofn = f'{odir}/{rx.split("/")[-1]}'
            download_request(url+rx, ofn, force=force, hr=hr)
    elif db == 'bev':
        rxlist = []
        if hr:
            return
        else:
            url = f'{urllist[db]}/obs/{year}/{doy}/'
            print (url)
            with urllib.request.urlopen(url) as response:
                html = response.read().decode('ascii')
                soup = BeautifulSoup(html, 'html.parser')
                for link in soup.find_all('a'):
                    if link.get('href') is not None and len(link.get('href')) == 15 and link.get('href').endswith("d.gz"):
                        rxlist.append(link.get('href')) 
                    elif link.get('href') is not None and link.get('href').endswith("MO.crx.gz"):
                        rxlist.append(link.get('href')) 
            rxlist = np.array(rxlist)
            # rxpath = np.array(rxpath)
            rxnames = np.array([r[:4].lower() for r in rxlist])
            if isinstance(rx, str):
                irx = np.isin(rxnames, rx)
                rxlist = rxlist[irx] if np.sum(irx) > 0 else None
                rxpath = rxpath[irx] if np.sum(irx) > 0 else None
            
            #    Download the data
            print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
            for rx in rxlist:
                ofn = f'{odir}{os.sep}{rx}'
                download_request(url+rx, ofn, force=force, hr=hr)
                
    elif db == 'brasil':
        if hr:
            print ("CORS does not support highrate data files")
            return 
        url = f'{urllist[db]}/{year}/{doy}/'
        rxlist = []
        r = requests.get(url, verify=False)
        if r.status_code == requests.codes.ok:
            for data in r:
                soup = BeautifulSoup(data.decode('ascii', 'ignore'), 'html.parser')
                for link in soup.find_all('a'):
                    if link.get('href') is not None and len(link.get('href')) == 12:
                        rxlist.append(link.get('href')[:4])
                    
        if isinstance(rx, str):
            irx = np.isin(np.asarray(rxlist), rx)
            rxlist = list(np.asarray(rxlist)[irx]) if np.sum(irx) > 0 else None
        
            
        if rxlist is not None:
            print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
            for rx in rxlist:
                path = f"{url}/{rx}{doy}1.zip"
                ofn = f'{odir}{rx}{doy}1.zip'
                download_request(urlpath=path, filename=ofn, force=force, hr=hr)
        else:
            print ('{} wasnt found'.format(rx))
            
    elif db in ('nz', 'nl', 'd'):
        if hr:
            print ("New Zealand DB doesnt have highrate data")
            return
        url = f'{urllist[db]}/{year}/{doy}/'
        rxlist = []
        r = requests.get(url, verify=False)
        if r.status_code == requests.codes.ok:
            for data in r:
                soup = BeautifulSoup(data.decode('ascii', 'ignore'), 'html.parser')
                for link in soup.find_all('a'):
                    if link.get('href') is not None and (link.get('href').endswith(('d.Z')) or link.get('href').endswith(('d.gz')) or link.get('href').endswith(('MO.crx.gz'))):
                        rxlist.append(link.get('href'))
                    
        rxlist = np.array(rxlist)
        # rxpath = np.array(rxpath)
        rxnames = np.array([r[:4].lower() for r in rxlist])
        if isinstance(rx, str):
            irx = np.isin(rxnames, rx)
            rxlist = rxlist[irx] if np.sum(irx) > 0 else None
            rxpath = rxpath[irx] if np.sum(irx) > 0 else None
        
        #    Download the data
        print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
        for rx in rxlist:
            ofn = f'{odir}{os.sep}{rx}'
            download_request(url+rx, ofn, force=force, hr=hr)
            
    elif db == 'es':
        if hr:
            print ("New Zealand DB doesnt have highrate data")
            return
        url = f'{urllist[db]}/{year}/{dt.strftime("%Y%m%d")}/'
        rxlist = []
        r = requests.get(url, verify=False)
        if r.status_code == requests.codes.ok:
            for data in r:
                soup = BeautifulSoup(data.decode('ascii', 'ignore'), 'html.parser')
                for link in soup.find_all('a'):
                    if link.get('href') is not None and (link.get('href').endswith(('d.Z')) or link.get('href').endswith(('MO.crx.gz'))):
                        rxlist.append(link.get('href'))
                    
        rxlist = np.array(rxlist)
        # rxpath = np.array(rxpath)
        rxnames = np.array([r[:4].lower() for r in rxlist])
        if isinstance(rx, str):
            irx = np.isin(rxnames, rx)
            rxlist = rxlist[irx] if np.sum(irx) > 0 else None
            rxpath = rxpath[irx] if np.sum(irx) > 0 else None
        
        #    Download the data
        print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
        for rx in rxlist:
            ofn = f'{odir}{os.sep}{rx}'
            download_request(url+rx, ofn, force=force, hr=hr)
    elif db == 'epos':
        if hr:
            print ("New Zealand DB doesnt have highrate data")
            return
        if dt.year < 2016:
            url = f'{urllist[db]}/RINEX2/30s/{year}/{doy}/'
        else:
            url = f'{urllist[db]}/RINEX3/30s/{year}/{doy}/'
        
        rxlist = []
        r = requests.get(url, verify=False)
        if r.status_code == requests.codes.ok:
            for data in r:
                soup = BeautifulSoup(data.decode('ascii', 'ignore'), 'html.parser')
                for link in soup.find_all('a'):
                    if link.get('href') is not None and (link.get('href').endswith(('D.Z')) or link.get('href').endswith(('MO.crx.gz'))):
                        rxlist.append(link.get('href'))
                    
        rxlist = np.array(rxlist)
        rxnames = np.array([r[:4].lower() for r in rxlist])
        if isinstance(rx, str):
            irx = np.isin(rxnames, rx)
            rxlist = rxlist[irx] if np.sum(irx) > 0 else None
            rxpath = rxpath[irx] if np.sum(irx) > 0 else None
        
        #    Download the data
        print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
        for rx in rxlist:
            ofn = f'{odir}{os.sep}{rx}'
            download_request(url+rx, ofn, force=force, hr=hr)
            
    elif db in ('sonel', 'fr', 'fr2'):
        ftp = ftplib.FTP(urllist[db])
        ftp.login()
        if hr:
            print ("Sonel database does not support highrate data quite yet")
            return
        else:
            if db == 'sonel':
                rpath = f'gps/data/{year}/{doy}/'
            elif db == 'fr':
                rpath = f'pub/data_v3/{year}/{doy}/data_30/'
            elif db == 'fr2':
                rpath = f'data/{year}/{doy}'
        ftp.cwd(rpath)
        rxlist = np.array(getStateList(year, doy, ftp, db, rxn=rx))
        if isinstance(rx, str):
            irx = np.isin(np.asarray(rxlist), rx)
            rxlist = list(np.asarray(rxlist)[irx]) if np.sum(irx) > 0 else None
        
        print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
        for urlrx in rxlist:
            download_cddis(ftp, urlrx, odir+urlrx, force=force)
    
    elif db == 's':
        ftp = ftplib.FTP(urllist[db])
        ftp.login('smrak', 'e_8m1]t;#%')
        if hr:
            print ("Sonel database does not support highrate data quite yet")
            return
        rpath = 'rinex3/{year}/{doy}/'
        rxlist = np.array(getStateList(year, doy, ftp, db, rxn=rx))
        if isinstance(rx, str):
            irx = np.isin(np.asarray(rxlist), rx)
            rxlist = list(np.asarray(rxlist)[irx]) if np.sum(irx) > 0 else None
        
        print ('Downloading {} receivers to: {}'.format(len(rxlist), odir))
        for urlrx in rxlist:
            download_cddis(ftp, urlrx, odir+urlrx, force=force)
            
    elif db == 'ring':
        if hr:
            print ("CORS does not support highrate data files")
            return 
        url = urllist[db]
        ftp = ftplib.FTP(url)
        ftp.login()
        rpath = f'/OUTGOING/RINEX30/RING/{year}/{doy}/'
        ftp.cwd(rpath)
        rxlist = np.array(getStateList(year, doy, ftp, db, rxn=rx))
        if isinstance(rx, str):
            irx = np.isin(np.asarray(rxlist), rx)
            rxlist = list(np.asarray(rxlist)[irx]) if np.sum(irx) > 0 else None
        for urlrx in rxlist:
            download_cddis(ftp, urlrx, odir+urlrx, force=force)
    else:
        raise('Wrong database')
if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('date', help='2017-5-27, or 2017-251, or a comma-separed start,end dates', type=str)
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
    dates = P.date.split(',')
    if len(dates) < 1:
        print("the date is a required argument")
        exit()
    elif len(dates) == 2:
        dtdates = np.arange(parser.parse(dates[0]),parser.parse(dates[1])+timedelta(hours=1), timedelta(days=1)).astype('datetime64[s]').astype(datetime)
        dates = [t.strftime("%Y-%m-%d") for t in dtdates]
    elif len(dates) > 2:
        print("Can't compile the dates")
        exit()
    for d in dates:
        if P.db == 'all':
            a = ['cors', 'cddis', 'unavco', 'brasil', 'chain', 'euref', 'ring', 'bev', 'sonel', 'nz', 'nl', 'd', 'fr', 'epos', 'es']
            for db in a:
                try:
                    getRinexObs(date = d, db = db, 
                            odir = P.dir, rx = P.rx, dllist = P.dllist, 
                            hr = P.highrate, force = P.force, fix = P.fixpath)
                except:
                    print (f"{db} Didn't work")
        elif P.db == 'highlat':
            a = ['cors', 'cddis', 'unavco', 'euref', 'sonel', 'epos']
            for db in a:
                getRinexObs(date = d, db = db, 
                            odir = P.dir, rx = P.rx, dllist = P.dllist, 
                            hr = P.highrate, force = P.force, fix = P.fixpath)
        elif P.db == 'conus':
            a = ['cors', 'unavco', 'cddis']
            for db in a:
                getRinexObs(date = d, db = db, 
                            odir = P.dir, rx = P.rx, dllist = P.dllist, 
                            hr = P.highrate, force = P.force, fix = P.fixpath)
        elif P.db == 'south':
            a = ['cddis', 'unavco', 'brasil', 'sonel']
            for db in a:
                getRinexObs(date = d, db = db, 
                            odir = P.dir, rx = P.rx, dllist = P.dllist, 
                            hr = P.highrate, force = P.force, fix = P.fixpath)
        elif P.db == 'europe':
            a = ['euref', 'ring', 'bev', 'sonel', 'epos', 'd', 'nl', 'es', 'fr']
            for db in a:
                getRinexObs(date = d, db = db, 
                            odir = P.dir, rx = P.rx, dllist = P.dllist, 
                            hr = P.highrate, force = P.force, fix = P.fixpath)
        else:
            getRinexObs(date = d, db = P.db, odir = P.dir, 
                        rx = P.rx, dllist = P.dllist, hr = P.highrate, 
                        force = P.force, fix = P.fixpath)
            
