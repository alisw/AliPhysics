#!/bin/sh

# 1 argument      - start file
# 2 argument      - number of files
# 3 argument      - input file list
# 4 argument      - run number used for ConfigOCDB

#myvar=0;
#while [ $myvar -ne 640 ] ; do bsub -q alice-t3_8h  submitcalib.sh $myvar $(($myvar+1))  /lustre_alpha/alice/miranov/rec/october2008/all/esdzip.txt; myvar=$(( $myvar +2 )) ; echo $myvar ; done

echo Hallo world
echo Hostname $HOSTNAME
df /tmp

echo $1   $2  $3
source ../balice.sh
source ../alienSetup.sh

mkdir $1_$2
cp *.C $1_$2 
cd $1_$2
cp $3 esd.txt

fstart=$1
fend=$2
cdbrun=$4
echo start aliroot
echo command aliroot -q -b CalibrateTPC.C\($fstart,$fend,$cdbrun\)
echo PWD `pwd`
command aliroot  -q -b "CalibrateTPC.C($fstart,$fend,$cdbrun)"
rm list.txt


