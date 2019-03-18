#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 11:55:08 2019

@author: smrak
"""

import numpy as np
import georinex as gr
import os

# Files
NCF = '/media/smrak/gnss/obs/2017/147/ab011470.17d.nc'
NF  = '/media/smrak/gnss/nav/brdc1470.17n'
SBF = '/media/smrak/gnss/jplg/jplg1470.17i'

SB = True
el_mask = 30
tlim=None
H = 350


svlist = gr.load(NCF).sv.values
navdata = gr.load(NF)
navdatatime = navdata.time.values