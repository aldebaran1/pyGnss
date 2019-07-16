#!/usr/bin/env bash

NAVDIR=/media/smrak/gnss/nav/
SBDIR=/media/smrak/gnss/jplg/
OBSDIR=/media/smrak/gnss/obs/

if [ $1 ]; then
  date=$1
else
  echo Enter the date YYYY-mm-dd:
  read date
fi

ddt=$(date -d $date +'%Y-%j')
year=$(date -d $date +'%Y')
doy=$(date -d $date +'%j')
foldername=$(date -d $date +'%m%d')

FN=$OBSDIR$year/$foldername/
FNYAML=$OBSDIR$year/all$foldername.yaml
FNYAMLC=$OBSDIR$year/conus$foldername.yaml

echo Downloading 1/7:
python download_rnxi.py $ddt $SBDIR
python download_rnxn.py $ddt $NAVDIR
python download_rnxn.py $ddt $NAVDIR --type sp3
python download_rnxo.py $ddt conus $OBSDIR
echo Unzipping 2/7:
python unzip_rm.py $FN
echo ==============================
echo listing receiver locations 3/7 from $FN
python rxlist.py $FN
echo ==============================
echo filter receiver locations 4/7
echo $FN'rxlist'$foldername.$year.h5 $FNYAML
python rxfilter.py $FN'rxlist'$foldername.$year.h5 $FNYAML
python rxfilter.py $FN'rxlist'$foldername.$year.h5 $FNYAMLC --lonlim -140 -50 --latlim 0 65
echo ==============================
echo Converting 5/7
python rnx2nc.py $FN
echo Deleting unconverted and rnx/crx
rm -rf $FN*.*d
rm -rf $FN*.*o
