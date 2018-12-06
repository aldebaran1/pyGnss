# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 10:50:06 2018

@author: smrak
"""

import os
from glob import glob

def getSatBias(fn: str = None):
    """
    Read satellite biases in TECu from JPL .XXi rinex file
    ftp://cddis.nasa.gov/gnss/products/ionex/
    Returns a dictionaty svbias, {'G05': -8.2} 
    """
    if os.path.isdir(fn):
        fn = glob(fn + '*.**i')
        if len(fn) >= 1:
            fn = fn[0] 
        else: 
            raise('No satbias files found in the folder')
            return 0
    i = 0
    svbias = {}
    with open(fn, 'r') as f:
        while 1:
            line = f.readline()
            if 'DIFFERENTIAL CODE BIASES' in line:
                while i < 32:
                    l = f.readline().split()[:2]
                    svbias['G'+l[0]] = float(l[1])
                    i += 1
                f.close()
                break
    return svbias