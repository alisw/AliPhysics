#!/bin/bash

# Script to run:
#    1. reconstruction
#    2. calibration and friend track filtering
#    3. tag creation
#
# Files assumed to be in working directory:
# rec.C               - reconstruction macro
# runCalibTrain.C     - calibration/filtering macro
# Arguments:
#    1  - raw data file name
#    2  - number of events to be processed
#    3  - run number 

# example:
# runPassX.sh raw.root  50  104892

#ALIEN setting
entries=1000
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
echo ">>>>>>>>> rec.C is..."
cat rec.C
echo

echo ">>>>>>> Running AliRoot to reconstruct $1. Run number is $runnum..."
aliroot -l -b -q rec.C\(\"alien://$1\"\) 2>&1 | tee rec.log

echo ">>>>>>> Running AliRoot to make calibration..."
aliroot -l -b -q  runCalibTrain.C\(\"$runnum\"\)   2>&1 | tee calib.log

echo ">>>>>>> Running AliRoot to generate Tags..."
aliroot -l -b -q tag.C\(\) 2>&1 | tee tag.log
