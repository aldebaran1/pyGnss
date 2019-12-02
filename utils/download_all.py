#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 12:34:41 2019

@author: smrak
"""

import subprocess
import multiprocessing as mp
from dateutil import parser
from datetime import timedelta
from time import sleep
import numpy as np
from argparse import ArgumentParser

def _getDates(start, stop):
    d = parser.parse(start)
    dates = [start]
    while d < parser.parse(stop):
        d += timedelta(days=1) 
        dates.append(d.strftime("%Y-%m-%d"))
    
    return np.array(dates)

def startP(line):
    subprocess.call(line, shell=True)
    sleep(1)
    return "Done with {}".format(line)
    
def main(start, stop, root=None, p=1):
    if root is None:
        root = './batch_tid.sh'
    dates_list = _getDates(start, stop)
    cmmnds = [root + ' ' + d for d in dates_list]
    
    pool = mp.Pool(processes=p)
    results = [pool.apply_async(startP, args=(l,)) for l in cmmnds]
    output = [p.get() for p in results]
    print (output)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('start', type = str, help='Start date "yyyy-mm-dd"')
    p.add_argument('stop', type = str, help='Stop date "yyyy-mm-dd"')
    p.add_argument('--program', default=None, type=str, help='Which bash program to run? Defult=batch_tid.sh')
    p.add_argument('-j', '--core', type=int, help='Define number of parllel processes', default=1)
    P = p.parse_args()
    
    main(P.start, P.stop, root=P.program, p=P.core)
    
