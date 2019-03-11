# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 13:34:38 2018

@author: smrak
"""

from glob import glob

folder = 'C:\\Users\\smrak\\Google Drive\\BU\\projects\\Themahali\\data\\rinex\\2015-10-07\\'

fn = glob(folder + '*.**i')[0]

i = 0
svbias = []
with open(fn, 'r') as f:
    while 1:
        line = f.readline()
        if 'DIFFERENTIAL CODE BIASES' in line:
            while i < 32:
                l = f.readline().split()[:2]
                b = ('G'+l[0], float(l[1]))
                svbias.append(b)
                i += 1
            f.close()
            break