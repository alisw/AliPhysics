#!/bin/bash

# Script to run Calibration Summary Extraction for the TPC :
# Arguments:
#    1  - run number 

# example:
# runCalibSummary.sh 119037

runnum=$1

echo
echo Run number to be processed $runnum
echo
echo ===========================
echo ">>>>>>>>> PATH is..."
echo $PATH
echo ">>>>>>>>> ROOTSYS is..."
echo $ROOTSYS
echo ">>>>>>>>> LD_LIBRARY_PATH is..."
echo $LD_LIBRARY_PATH
echo ==========================
echo

echo ">>>>>>> Running AliRoot to extract calibration summary..."
aliroot -l -b -q  ./runCalibSummary.C\(\"$runnum\"\)   2>&1 | tee calib.log

if [ -f dcsTime.root ]
    then
    mv dcsTime.root calibSummary_run$runnum.root
else
    echo
    echo "dcsTime.root file was not created - job failure"
fi
