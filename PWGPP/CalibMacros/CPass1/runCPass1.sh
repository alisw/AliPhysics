#!/bin/bash

# Script to run:
#    1. reconstruction
#    2. calibration and friend track filtering
#
# Files assumed to be in working directory:
# recCPass1.C          - reconstruction macro
# runCalibTrain.C     - calibration/filtering macro
# Arguments (run locally):
#    1  - raw data file name
#    2  - number of events to be processed
#    3  - run number 

# example:
# runCPass1.sh raw.root  50  104892

#ALIEN setting
# $1 = raw input filename
runNum=`echo $1 | cut -d "/" -f 6 | sed 's/^0*//'`
if [ $# -eq 1 ] ; then
  # alien Setup
  nEvents=99999999
  fileName="alien://"$1
  ocdbPath="raw://"
fi;
if [ $# -eq 4 ] ; then
  # local setup
  nEvents=$2
  runNum=$3
  fileName=$1
  ocdbPath=$4
fi

echo xxxxxxxxxxxxxxxxxxxxxxxxxxx
echo runCPass1.sh Input arguments
echo fileName=$fileName
echo nEvents=$nEvents
echo runNum=$runNum
echo ocdbPath=$ocdbPath
echo xxxxxxxxxxxxxxxxxxxxxxxxxxx

if [ -f Run0_999999999_v3_s0.root ]; then
    mkdir -p TPC/Calib/Correction
    mv Run0_999999999_v3_s0.root TPC/Calib/Correction/
fi



echo File to be  processed $1
echo Number of events to be processed $nEvents

echo ">>>>>>>>> PATH is..."
echo $PATH
echo ">>>>>>>>> LD_LIBRARY_PATH is..."
echo $LD_LIBRARY_PATH
echo ">>>>>>>>> recCPass1.C is..."
#cat recCPass1.C
echo

echo ">>>>>>> Running AliRoot to reconstruct $1. Run number is $runNum..."

mkdir Barrel
mkdir OuterDet

cd Barrel
time aliroot -l -b -q ../recCPass1.C\(\""$fileName\", $nEvents, \"$ocdbPath"\"\) 2>&1 | tee rec_Barrel.log
mv syswatch.log syswatch_rec_Barrel.log

echo ">>>>>>> Running AliRoot to make calibration..."
time aliroot -l -b -q ../runCalibTrain.C\(\""$runNum\",\"AliESDs.root\",\"$ocdbPath"\"\)   2>&1 | tee calib.log
mv syswatch.log syswatch_calib.log
echo ">>>>>>> Doing ls -l"
ls -l

echo ">>>>>>> Running the QA train..."
time aliroot -b -q ../QAtrain.C\($runNum\) 2>&1 | tee qa_Barrel.log

for file in *.stat; do
    mv $file $file.qa
done

echo ">>>>>>>> Moving files to upper directory"
mv AliESDs.root ../AliESDs_Barrel.root
mv rec_Barrel.log ../rec_Barrel.log
mv calib.log ../calib.log
mv AliESDfriends_v1.root ../AliESDfriends_v1.root
mv qa_Barrel.log ../qa_Barrel.out
mv QAresults.root ../QAresults_Barrel.root
if [ -f AODtpITS.root ] ; then
 mv AODtpITS.root ../
fi

cd ../OuterDet
time aliroot -l -b -q ../recCPass1_OuterDet.C\(\""$fileName\", $nEvents, \"$ocdbPath"\"\) 2>&1 | tee rec_Outer.log
mv syswatch.log syswatch_rec_Outer.log

echo ">>>>>>> Running the QA train..."
time aliroot -b -q ../QAtrain.C\($runNum\) 2>&1 | tee qa_Outer.log

for file in *.stat; do
    mv $file $file.qa
done

mv AliESDs.root ../AliESDs_Outer.root
mv rec_Outer.log ../rec_Outer.log
mv qa_Outer.log ../qa_Outer.out
mv QAresults.root ../QAresults_Outer.root
