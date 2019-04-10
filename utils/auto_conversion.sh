#!/usr/bin/env bash
X0=-140
X1=-55
Y0=10
Y1=57
H=350
CFG=/home/smrak/Documents/pyTID/pytid/cfg/t2g_conus.yaml

NAVDIR=/media/smrak/gnss/nav/
SBDIR=/media/smrak/gnss/jplg/
OBSDIR=/media/smrak/gnss/obs/

PYCONV=/home/smrak/Documents/pyTID/pytid/singlerx_v2.py
PYGRID=/home/smrak/Documents/pyTID/pytid/tid2grid_v2.py

echo Enter the date YYYY-mm-dd:
read date
ddt=$(date -d $date +'%Y-%j')
year=$(date -d $date +'%Y')
doy=$(date -d $date +'%j')

FN=/media/smrak/gnss/obs/$year/$doy/
FNYAML=/media/smrak/gnss/obs/$year/conus$doy.yaml
FNGRID=/media/smrak/gnss/hdfgrid/
HDFFN=conus$doy.yaml_$year'_'$H'km.h5'
echo Downloading 1/7:
python download_rnxi.py $ddt $SBDIR
python download_rnxn.py $ddt $NAVDIR
python download_rnxo.py $ddt conus $OBSDIR
echo Unzipping 2/7:
python unzip_rm.py $FN
echo ==============================
echo listing receiver locations 3/7 from $FN
python rxlist.py $FN
echo ==============================
echo filter receiver locations 4/7
echo $FN'rxlist'$doy.$year.h5 $FNYAML $X0 $X1 $Y0 $Y1
python rxfilter.py $FN'rxlist'$doy.$year.h5 $FNYAML $X0 $X1 $Y0 $Y1
echo ==============================
echo Converting 5/7
python rnx2nc.py $FN --list $FNYAML
echo ==============================
echo Processing 6/7
python $PYCONV $year $doy $FNYAML --altkm $H
echo ==============================
echo Gridding 7/7
python $PYGRID $HDFFN $FNGRID'conus_03_'$year'_'$doy.h5
