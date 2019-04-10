#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 10:38:16 2018

@author: Sebastijan Mrak <smrak@bu.edu>
"""

import glob
import os
import subprocess
import platform
from numpy import unique

def unzip(f, delete=False, timeout=30):
    if platform.system() == 'Linux':
        try:
            print('Decompressing: ', f)
            subprocess.call('./CRX2RNX ' + f, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)
            if delete:
                print('Deleting: ', f)
                subprocess.call('rm -r ' + f, shell=True)
        except:
            print('Problems with: ',f)
    elif platform.system() == 'Windows':
        try:
            print('Decompressing: ', f)
            subprocess.call('crx2rnx.exe "{}"'.format(f), shell=True, timeout=timeout)
            if delete:
                print('Deleting: ',f)
                subprocess.call('del "{}"'.format(f), shell=True,timeout=timeout)
        except:
            print('Problems with: ',f)
        
    return

def crx2rx(folder, delete=False, timeout=10):
    """
    A function utilizin CRX2RNX script to decompress the RINEX observation 
    files with a *.YYd extension into a *.YYo extension. 
    The script also removes all unnecessary junk from the given directory e.g.,
    .crx and .YYd files.
    The input is a directory, and it automatically finds all files needed to be
    decompressed.
    """
    suffix = ['*.**d']
    for wlstr in suffix:
        if os.path.isdir(folder):
            filestr = os.path.join(folder,wlstr)
            flist = unique(sorted(glob.glob(filestr)))
            print ('Decompressing {} files in the folder'.format(len(flist)))
            for f in flist:
                unzip(f, delete=delete, timeout=timeout)
        else:
            unzip(folder, delete=delete, timeout=timeout)
    
    if delete:
        suffix.append('*.crx')
        for suffix in suffix:
            filestr = os.path.join(folder,suffix)
            flist = sorted(glob.glob(filestr))
            for f in flist:
                if platform.system() == 'Linux':
                    print('Deleting: ', f)
                    subprocess.call('rm -r ' + f, shell=True)
                elif platform.system() == 'Windows':
                    print('Deleting: ', f)
                    subprocess.call('del "{}"'.format(f), shell=True,timeout=15)

        
if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('folder',type=str)
    p.add_argument('-d', '--delete', help='Delete the resitudal .D file', action='store_true')
    p.add_argument('-t', '--tout', help='Timeout for decompressing', default=10, type=int)
    P = p.parse_args()
    
    crx2rx(P.folder, delete=P.delete, timeout=P.tout)