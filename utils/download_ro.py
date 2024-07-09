#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 16:53:38 2024

@author: mraks1
"""

import urllib.request
from bs4 import BeautifulSoup
from dateutil import parser
import numpy as np
import ssl, sys, subprocess, platform, wget, tarfile, os
from datetime import datetime, timedelta

def bar_progress(current, total, width=80):
  progress_message = "Downloading: %d%% [%d / %d] bytes" % (current / total * 100, current, total)
  # Don't use print() as it will print in new line every time.
  sys.stdout.write("\r" + progress_message)
  sys.stdout.flush()

def main(date, odir, const=None):
    # This restores the same behavior as before.
    context = ssl._create_unverified_context()
    ssl._create_default_https_context = ssl._create_unverified_context
    date = date.split(",")
    if len(date) == 1:
        dt = parser.parse(date)
        dates = [dt]
    elif len(date) == 2:
        dt = [parser.parse(d) for d in date]
        dates = np.arange(dt[0], dt[1]+timedelta(days=1), timedelta(days=1)).astype('datetime64[s]').astype(datetime)
    else:
        raise("Something went wrong with provided dates")
    for date in dates:
        year = date.year
        doy = date.strftime("%j")
        if const == 'cosmic2':
            url = f'https://data.cosmic.ucar.edu/gnss-ro/cosmic2/nrt/level1b/{year}/{doy}/'
            pattern = f'podTc2_nrt_{year}_{doy}.tar.gz'
        
            datalist = []
            with urllib.request.urlopen(url, context=context) as response:
                html = response.read().decode('ascii')
                soup = BeautifulSoup(html, 'html.parser')
                for link in soup.find_all('a'):
                    # if len(link.get('href')[:-1]) == 6:
                    datalist.append(link.get('href'))
                if pattern in np.array(datalist):
                    if not os.path.exists(odir):
                        if platform.system() == 'Windows':
                            subprocess.call(f'mkdir "{odir}"', shell=True)
                        else:
                            subprocess.call(f'mkdir -p "{odir}"', shell=True)
                    print (f"Downloading {url}{pattern}")
                    wget.download(url+pattern, out=odir, bar=bar_progress)
                    print (f"Unzipping {odir}{pattern}")
                    tar = tarfile.open(odir+pattern, "r:gz")
                    tar.extractall()
                    tar.close()
    
    return 0

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('date', help='Single date or start,end dates, comma separated')
    p.add_argument('odir', type=str, help='Oudput Directory')
    p.add_argument('-c', type=str, help="Which constallation? Supported: cosmic2", default='cosmic2')
    
    P = p.parse_args()
    
    main(P.date, P.odir)
    
    main(P.date, odir=P.odir)