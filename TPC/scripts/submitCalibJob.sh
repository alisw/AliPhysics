#!/bin/sh

# 1 argument      - start file
# 2 argument      - number of files
# 3 argument      - input file list
# 4 argument      - run number used for ConfigOCDB

fstart=$1
fend=$2
cdbrun=$4
run=$4
echo $1   $2  $3 $4

echo Hallo world
echo Hostname $HOSTNAME
echo RUN $4
df /tmp
source ../balice.sh
source ../alienSetup.sh


mkdir $1_$2
cp *.C $1_$2 
cd $1_$2
cp $3 esd.txt
#
#
#
dname=`basename \`pwd\``
basename=`pwd`/
tmpname=/tmp/$USER/$run/$dname
mkdirhier $tmpname
rm -rf  $tmpname/*
cp * $tmpname/
ls -al $tmpname
cd $tmpname
echo Working directory  $tmpname


echo start aliroot
echo command aliroot -q -b CalibrateTPC.C\($fstart,$fend,$cdbrun\)
echo PWD `pwd`
echo 00000000000000000000000000000000000000000000000
echo Start job
echo Date  `date`
echo 00000000000000000000000000000000000000000000000
command aliroot  -q -b "CalibrateTPC.C($fstart,$fend,$cdbrun)"
echo 00000000000000000000000000000000000000000000000
echo Copy output
echo Date  `date`
echo 00000000000000000000000000000000000000000000000
cp -rf * $basename
find $basename/ | grep .root

echo 00000000000000000000000000000000000000000000000
echo DELET OUTPY
echo Date  `date`
echo 00000000000000000000000000000000000000000000000
rm -rf  $tmpname/*
rm *.root


echo 00000000000000000000000000000000000000000000000
echo End job
echo Date `date`
echo 00000000000000000000000000000000000000000000000
exit


