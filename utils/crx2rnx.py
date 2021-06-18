#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 10:38:16 2018

@author: Sebastijan Mrak <smrak@bu.edu>
"""

import shutil
from pathlib import Path

import hatanaka


def crx2rx(folder, delete=False):
    """
    A function utilizing CRX2RNX script to decompress the RINEX observation
    files with a *.YYd extension into a *.YYo extension.
    The script also removes all unnecessary junk from the given directory e.g.,
    .crx and .YYd files.
    The input is a directory, and it automatically finds all files needed to be
    decompressed.
    """
    folder = Path(folder)
    suffix = ['*.??d']
    for wlstr in suffix:
        if folder.is_dir():
            flist = sorted(set(folder.glob(wlstr)))
            print(f'Decompressing {len(flist)} files in the folder')
            for f in flist:
                hatanaka.decompress_on_disk(f, delete=delete)
        else:
            hatanaka.decompress_on_disk(folder, delete=delete)

    if delete:
        suffix.append('*.crx')
        for suffix in suffix:
            for f in sorted(folder.glob(suffix)):
                print('Deleting: ', f)
                shutil.rmtree(f)


if __name__ == '__main__':
    from argparse import ArgumentParser

    p = ArgumentParser()
    p.add_argument('folder', type=str)
    p.add_argument('-d', '--delete', help='Delete the residual .D file', action='store_true')
    P = p.parse_args()

    crx2rx(P.folder, delete=P.delete)
