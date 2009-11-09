#!/bin/sh

# 1 argument      - start file
# 2 argument      - number of files
# 3 argument      - input file list
# 4 argument      - run number used for ConfigOCDB

echo Hallo world
echo Hostname $HOSTNAME
df /tmp

echo $1   $2  $3
source ../balice.sh
source ../alienSetup.sh
moval
mkdir $1_$2
cp *.C $1_$2 
cd $1_$2
cp $3 esd.txt

dname=`basename \`pwd\``
basename=`pwd`/${fstart}_${fend}
mkdirhier $basename
tmpname=/tmp/$USER/$dname/${fstart}_${fend}
mkdirhier $tmpname
cp * $tmpname
ls -al $tmpname
cd $tmpname



fstart=$1
fend=$2
cdbrun=$4
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

rm *.root
cp -rf * $basename


echo 00000000000000000000000000000000000000000000000
echo End job
echo Date `date`
echo 00000000000000000000000000000000000000000000000
exit


