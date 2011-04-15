#!/bin/bash

# Script to merge objects coming out of the calibration train:
# Arguments:
#    1 - directory on which to look for the files to be merged 

# example:
# mergeCalibObjects.sh /alice/cern.ch/user/z/zampolli/CalibTrain/output/

#ALIEN setting
# $1 = directory where to perform the find 
#runnum=$1

#echo Directory to look into = $1

echo ">>>>>>> Running AliRoot to merge calib objects found in $1"
aliroot -l -b -q merge.C\(\"$1\"\) 2>&1 | tee merge.log


