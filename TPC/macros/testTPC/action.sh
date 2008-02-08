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


mkdir -p /tmp/$USER/$1
cd /tmp/$USER/$1
which  xrdcp

command xrdcp $2 data.root

echo aliroot -b -q $ALICE_ROOT/TPC/macros/testTPC/$4.C

command aliroot -b -q $ALICE_ROOT/TPC/macros/testTPC/$4.C

for name in `ls *.root`; do
  xrdcp $name $3/$4/$name
done



