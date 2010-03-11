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
#entries=1000
# $1 = raw input filename
#runnum=`echo $1 | cut -d "/" -f 6`

#Local setting : setting variables 

entries=$2
runnum=$3
source $HOME/alienSetup.sh

echo File to be  processed $1
echo Number of events to be processed $entries
echo Run mumber $runnum

echo ALICE_ROOT = $ALICE_ROOT
echo AliROOT = $AliROOT
cp $ALICE_ROOT/.rootrc ~/.rootrc
cp $ALICE_ROOT/.rootrc $HOME
#cat $HOME/.rootrc
export GRID_TOKEN=OK

echo ">>>>>>>>> PATH is..."
echo $PATH
echo ">>>>>>>>> LD_LIBRARY_PATH is..."
echo $LD_LIBRARY_PATH
echo ">>>>>>>>> rec.C is..."
cat rec.C
echo


echo
echo ">>>>>>> Running AliRoot to reconstruct $1. Run number is $runnum..."
echo
if [ -e AliESDs.root ]; then
    echo AliESDs.root exist
    ls -al AliESD*
else
    echo aliroot -l -b -q rec.C\(\"$1\",$2\) 2>&1 | tee rec.log
    aliroot -l -b -q rec.C\(\"$1\",$2\) 2>&1 | tee rec.log
    echo aliroot -l -b -q tag.C\(\) 2>&1 | tee tag.log
    aliroot -l -b -q tag.C\(\) 2>&1 | tee tag.log
fi

echo
echo ">>>>>>> Running AliRoot to make calibration..."
echo 
echo aliroot -l -b -q  runCalibTrain.C\($runnum\)   2>&1 | tee calib.log
aliroot -l -b -q  runCalibTrain.C\($runnum\)   2>&1 | tee calib.log

echo
echo ">>>>>>> Running AliRoot to generate Tags..."
echo
echo aliroot -l -b -q tag.C\(\) 2>&1 | tee tag.log
#aliroot -l -b -q tag.C\(\) 2>&1 | tee tag.log
