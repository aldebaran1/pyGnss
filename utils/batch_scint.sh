#!/bin/bash -l

#$ -l h_rt=48:00:00 

CONVDIR=/usr3/graduate/smrak/pyTID/pytid/

NAVDIR=/projectnb/semetergrp/scint/nav/
SBDIR=/projectnb/semetergrp/scint/jplg/
OBSDIR=/projectnb/semetergrp/scint/obs/
HDFDIR=/projectnb/semetergrp/scint/hdf/

if [ $1 ]; then
  date=$1
else
  echo Enter the date YYYY-mm-dd:
  read date
fi

ddt=$(date -d $date +'%Y-%j')
year=$(date -d $date +'%Y')
doy=$(date -d $date +'%j')
mmdd=$(date -d $date +'%m%d')

FN=$OBSDIR$year/$mmdd/
FNYAML=$OBSDIR$year/reduced$mmdd'_d1_r2'.yaml

echo Downloading 1/7:
python download_rnxi.py $ddt $SBDIR
python download_rnxn.py $ddt $NAVDIR
python download_rnxo.py $ddt unavco $OBSDIR --highrate
echo Unzipping 2/7:
python unzip_rm.py $FN
echo ==============================
echo listing receiver locations 3/7 from $FN
python rxlist.py $FN
echo ==============================
echo filter receiver locations 4/7
echo $FN'rxlist'$mmdd.$year.h5 $FNYAML
python rxfilter.py $FN'rxlist'$mmdd.$year.h5 $FNYAML --density 1 -r 2
echo ==============================
echo Converting 5/7
python rnx2nc.py $FN 
echo Deleting unconverted and rnx/crx
rm -rf $FN*.*d
rm -rf $FN*.*o
echo ==============================
echo Processing 6/7 to tid
python $CONVDIR'nc2tec_v1.py' $date $FNYAML --ts 1 --elmask 30 --log --stec --cfg $CONVDIR'scint_processing.yaml'
