#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 10:38:16 2018

@author: Sebastijan Mrak <smrak@bu.edu>
"""

import glob
import os
import subprocess
from datetime import datetime
import platform
import numpy as np

def do_teqc(folder, rx, doy, fout=None):
    """
    A function utilizin CRX2RNX script to decompress the RINEX observation 
    files with a *.YYd extension into a *.YYo extension. 
    The script also removes all unnecessary junk from the given directory e.g.,
    .crx and .YYd files.
    The input is a directory, and it automatically finds all files needed to be
    decompressed.
    """
    if fout is None:
        fout = folder
    if not os.path.exists(fout):
        subprocess.call(f"mkdir -p '{fout}'", shell=True)
    fxn = folder + f"{os.sep}{rx}{doy}*.*o"
    ext = os.path.splitext(glob.glob(fxn)[0])[1]
    #teqc = f"teqc"
    cmd = f"{folder}{os.sep}{rx}{doy}[a-w]*{ext}"
    cmda = f"{fout}{os.sep}{rx}{doy}0{ext}"
    subprocess.call(f"teqc {cmd} > {cmda}", shell=True)
    return
    
if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('folder',type=str)
    p.add_argument('--odir', help='Output directory? If None, default is input directory', default=None)
    P = p.parse_args()
    
    rxfiles = np.array(sorted(glob.glob(P.folder + '*.*o')))
    rxnames = np.array([os.path.split(f)[1][:4] for f in rxfiles])
    doys = np.unique(np.array([int(os.path.split(f)[1][4:7]) for f in rxfiles]).astype(int))
    _, idx = np.unique(rxnames, return_index=True)
    rxall = rxnames[idx]
    print (f"Receivers {np.unique(rxnames).size}: {list(np.unique(rxnames))}")
    print (f"Doys: {list(doys)}")
    for doy in doys:
        for rx in rxall:
            t0 = datetime.now()
            do_teqc(P.folder, rx, doy, fout=P.odir)
            print (f"{rx}, {doy} converted in {datetime.now()-t0}")
