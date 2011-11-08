#!/bin/bash

# Script to run:
#    1. reconstruction
#    2. calibration and friend track filtering
#
# Files assumed to be in working directory:
# recPass0.C          - reconstruction macro
# runCalibTrain.C     - calibration/filtering macro
# Arguments (run locally):
#    1  - raw data file name
#    2  - number of events to be processed
#    3  - run number 

# example:
# runPass0.sh raw.root  50  104892

#ALIEN setting
# $1 = raw input filename
runnum=`echo $1 | cut -d "/" -f 6`

#Local setting
#entries=$2
#runnum=$3

echo File to be  processed $1
echo Number of events to be processed $entries

echo ">>>>>>>>> PATH is..."
echo $PATH
echo ">>>>>>>>> LD_LIBRARY_PATH is..."
echo $LD_LIBRARY_PATH
echo ">>>>>>>>> recPass0.C is..."
cat recPass0.C
echo

echo ">>>>>>> Running AliRoot to reconstruct $1. Run number is $runnum..."
aliroot -l -b -q recPass0.C\(\"alien://$1\"\) 2>&1 | tee rec.log
mv syswatch.log syswatch_rec.log

echo ">>>>>>> Running AliRoot to make calibration..."
aliroot -l -b -q  runCalibTrain.C\(\"$runnum\"\)   2>&1 | tee calib.log
mv syswatch.log syswatch_calib.log

