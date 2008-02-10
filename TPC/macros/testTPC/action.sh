#!/bin/bash
#aguments
#1 TString jobID, 
#2 TString inputData
#3 TString outputDir
#4 TString   action

echo $1
echo $2
echo $3
echo $4


olddir=`pwd`

mkdir -p /tmp/$USER/$1
mkdir -p /tmp/$USER/$1/$4

cd /tmp/$USER/$1/$4
which  xrdcp

xrdcp -np $2 data.root

echo aliroot -b -q $ALICE_ROOT/TPC/macros/testTPC/$4.C

aliroot -b -q $ALICE_ROOT/TPC/macros/testTPC/$4.C
echo END ACTION $1

rm data.root
ls -al 
echo CREATING ZIP FILE
zip -n root $4 *.root
echo COPING DATA
for name in `ls *.root`; do
    echo  xrdcp $name $3/$4/$name
  xrdcp -np $name $3/$4/$name
done
echo xrdcp -np  $4.zip  $3/$4.zip
 xrdcp -np $4.zip  $3/$4.zip


cd $olddir

