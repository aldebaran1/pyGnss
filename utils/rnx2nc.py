#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 13:25:38 2017

@author: Sebastijan Mrak <smrak@gmail.com>
"""
import os
import glob
import platform
from numpy import array, isin, logical_not
import yaml
from time import sleep
import georinex as gr
from datetime import datetime

def _convert(file, odir, i, tlim=None, fast=True, use = 'G'):
    try:
        gr.load(file,out=odir, useindicators=i, 
                fast = fast, tlim=tlim, use = use)
    except Exception as e:
        print (e)
    sleep(0.1)
    
def _iterate(file, odir, override, indicators, tlim=None, use = 'G', fast = True):
    tail = os.path.split(file)[1]
    newfn = file + '.nc'
    if not override:
        if not os.path.exists(newfn):
            print ('Converting: {}'.format(tail))
            _convert(file, odir, indicators, tlim=tlim, fast=fast, use=use)
        else:
            print ("{} already exists.".format(os.path.split(newfn)[1]))
    else:
        print ('Converting: {}'.format(tail))
        _convert(file, odir, indicators, tlim=tlim, fast=fast, use=use)
        
def convertObs2HDF(folder=None, sufix=None, odir=None, override=False,
                   indicators=False, tlim=None, rxlist = None, use = 'G', fast = True):
    """
    This script converts RINEX 2.11 observation files in a given directory into
    a hdf5 organized data structure, utilizing pyRINEX script. Find the script
    in the main directory.
    """
    if odir is None:
        if os.path.isdir(folder):
            odir = folder
        elif os.path.isfile(folder):
            odir = os.path.split(folder)[0]
        else:
            raise('odir must be either a folder or a directory')
    if rxlist is not None:
        assert os.path.isfile(rxlist), 'Input argument must be a file, not folder!'
        assert os.path.splitext(rxlist)[1] == '.yaml', 'File must be .yaml format.'
        
        stream = yaml.load(open(rxlist, 'r'))
        rxl = array(sorted(stream.get('rx')))
        print('Number of files to convert: {}'.format(stream.get('total')))
        sx = '*.*d' if sufix is None else sufix
        flist = array(sorted(glob.glob(os.path.join(folder, sx))))
        fnamelist = array([os.path.split(r)[1][:4] for r in flist])
        idl = isin(fnamelist, rxl)
        t0 = datetime.now()
        print (override)
        for i, file in enumerate(flist[idl]):
            _iterate(file, odir, override=override, indicators=indicators, 
                     tlim=tlim, fast=fast, use=use)
        t1 = datetime.now()
        print ('Conversion successfull. It took {}s to convert them, at an average of {} per file.'.format(t1-t0, (t1-t0)/flist[idl].shape[0]))

    else:
        if os.path.isdir(folder):
            if sufix is None:
                wlist = ['*.**o', '*.**O', '*.**d']
            else:
                wlstr = sufix
            for wlstr in wlist:
                filestr = os.path.join(folder,wlstr)
                flist = sorted(glob.glob(filestr))
                for file in flist:
                    _iterate(file, odir, override, indicators, tlim=tlim, fast=fast, use=use)
        elif os.path.isfile(folder):
                file = folder
                _iterate(file, odir, override, indicators, tlim=tlim, fast=fast, use=use)
        else:
            raise ("Something went wrong, dude")
if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('folder', type = str)
    p.add_argument('-odir', '--odir', help = 'Destination folder, if None-> the same as input folder', default = None)
    p.add_argument('-f', '--force', help = "Force override, if the NC file already exist", action = 'store_true')
    p.add_argument('-i', '--indicators', help = "Parse & store the indicators (lli/ssi)?", action = 'store_true')
    p.add_argument('--fast', help = "Disable fast processing. Default == True", action = 'store_true')
    p.add_argument('-s', '--sufix', help = 'specify a sufix for desired observation files', type = str, default = None)
    p.add_argument('--use', help = 'Use specific constallations only. Default == "G"', type = str, default = 'G')
    p.add_argument('--tlim', help = 'set time limints for the file to cenvert', nargs = 2)
    p.add_argument('--list', help = 'Recuced list of receivers (.yaml dict) to convert to .NC', type = str, default = None)
    P = p.parse_args()
    
    convertObs2HDF(folder = P.folder, sufix=P.sufix, odir=P.odir, override=P.force, 
                   indicators=P.indicators, tlim=P.tlim, rxlist=P.list, fast=logical_not(P.fast), use=P.use)