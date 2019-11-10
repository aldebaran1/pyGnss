#!/usr/bin/env bash
X0=-180
X1=-50
Y0=0
Y1=90

NAVDIR=/media/smrak/gnss/nav/
SBDIR=/media/smrak/gnss/jplg/
OBSDIR=/media/smrak/gnss/obs/highrate/

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

FN=/media/smrak/gnss/obs/highrate/$year/$mmdd/
FNYAML=/media/smrak/gnss/obs/highrate/$year/reduced$mmdd'_d1_r2'.yaml

echo Downloading 1/5:
python download_rnxi.py $ddt $SBDIR
python download_rnxn.py $ddt $NAVDIR
python download_rnxo.py $ddt unavco $OBSDIR --highrate
echo Unzipping 2/5:
python unzip_rm.py $FN
echo ==============================
echo listing receiver locations 3/5 from $FN
python rxlist.py $FN
echo ==============================
echo filter receiver locations 4/5
echo $FN'rxlist'$mmdd.$year.h5 $FNYAML $X0 $X1 $Y0 $Y1
python rxfilter.py $FN'rxlist'$mmdd.$year.h5 $FNYAML --lonlim $X0 $X1 --latlim $Y0 $Y1 --density 1 -r 2
echo ==============================
echo Converting 5/5
python rnx2nc.py $FN --list $FNYAML
