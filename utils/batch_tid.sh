#!/usr/bin/env bash

CONVDIR=/usr3/graduate/smrak/pyTID/pytid/

NAVDIR=/projectnb/semetergrp/tid/nav/
SBDIR=/projectnb/semetergrp/tid/jplg/
OBSDIR=/projectnb/semetergrp/tid/obs/
HDFDIR=/projectnb/semetergrp/tid/hdf/

if [ $1 ]; then
  date=$1
else
  echo Enter the date YYYY-mm-dd:
  read date
fi

ddt=$(date -d $date +'%Y-%j')
year=$(date -d $date +'%Y')
doy=$(date -d $date +'%j')

FN=$OBSDIR$year/$doy/
FNYAML=$OBSDIR$year/all$doy.yaml
FNYAMLC=$OBSDIR$year/conus$doy.yaml

#echo Downloading 1/7:
#python download_rnxi.py $ddt $SBDIR
#python download_rnxn.py $ddt $NAVDIR
#python download_rnxo.py $ddt conus $OBSDIR
#echo Unzipping 2/7:
#python unzip_rm.py $FN
#echo ==============================
#echo listing receiver locations 3/7 from $FN
#python rxlist.py $FN
#echo ==============================
#echo filter receiver locations 4/7
#echo $FN'rxlist'$doy.$year.h5 $FNYAML
#python rxfilter.py $FN'rxlist'$doy.$year.h5 $FNYAML
#python rxfilter.py $FN'rxlist'$doy.$year.h5 $FNYAMLC --lonlim -140 -50 --latlim 0 65
#echo ==============================
#echo Converting 5/7
#python rnx2nc.py $FN 
#echo Deleting unconverted and rnx/crx
#rm -rf $FN*.*d
#rm -rf $FN*.*o
#echo ==============================
#echo Processing 6/7 to tid
echo python $CONVDIR'nc2tec_v1.py' $date $FNYAMLC --ts 120 --elmask 20 --log --cfg $CONVDIR'tid_processing.yaml'
