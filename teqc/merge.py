# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 15:07:59 2018

@author: smrak
"""
import subprocess
from glob import glob
import os
from argparse import ArgumentParser

def main(folder: str, rinexfile: str, ofn: str = None,
         verbose: bool = False):
#    folder = "D:\\RefData.17\\Month.Aug\\Day.21\\MOAT\\"
#    folder = "E:\\Day.21\\MOJC\\"
#    rinexfile = 'o'
    
    flist = sorted(glob(folder+'*'+rinexfile))
    
    head, fn = os.path.split(flist[0])
    outfn = fn[:7] + '0.'+fn[-3:-1]+rinexfile
    of = folder+outfn
    
    arg = " ".join(flist)
    
    teqc = folder + 'teqc.exe '
    #teqc = '"C:\\Users\\smrak\\Google Drive\\BU\\software\\pyGnss\\teqc\\teqc.exe" '
    teqcline = teqc + arg + ' > ' + of
    if verbose:
        print (teqcline)
    
    subprocess.call(teqcline, shell=True)
    

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('folder',type=str)
    p.add_argument('-t', '--type', help='Rinex type to convert? (o,n,g,..)', default='o')
    p.add_argument('-odir', '--odir', help='path to output director=?', default=None)
    p.add_argument('-v', '--verbose', help='verbose', action='store_true')
    P = p.parse_args()
    
    main(folder=P.folder, rinexfile=P.type, ofn=P.odir, verbose=P.verbose)