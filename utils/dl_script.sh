#!/bin/bash -l

HOME='/homes/mraks1/software/pyGnss/utils/'
ROOT='/disks/d0395/mraks1/rinex/'

#Load Python Env
#source /homes/mraks1/venv/bin/activate.csh

#Mount /d0395
cd /disks/d0395/
cd $HOME

if [ $1 ]; then
  date=$1
else
  echo Enter the date yyyy-mm-dd
  read date
fi

if [ $2 ]; then
  mode=$2
else
  mode='all'
fi
mode='all'

echo $date
echo $mode

# directory tree from date
year=$(date -d $date +'%Y')
mmdd=$(date -d $date +'%m%d')

ODIR=$ROOT$year/$mmdd/
NHP='nohup_dl_'$year$mmdd'.out'


# 1. Download rnxo
echo Downloading RNXO Data: 
python $HOME'download_rnxo.py' $date $mode $ODIR '--fixpath'
python $HOME'download_rnxn.py' $date $ROOT'nav/'
# 2. Unzip all the files
python $HOME'unzip_rm.py' $ODIR
# 3. List the files
python $HOME'rxlist.py' $ODIR
# 4. Filter the files
python $HOME'rxfilter.py' $ODIR'rxlist'$mmdd.$year.h5


