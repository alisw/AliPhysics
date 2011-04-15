#!/bin/bash

# Script to merge objects coming out of the calibration train:
# Arguments:
#    1 - directory on which to look for the files to be merged 
#    2 - run number
#    3 - OCDB output path

# example:
# mergeMakeOCDB.sh /alice/cern.ch/user/z/zampolli/CalibTrain/output/ 000114798 alien://folder=/alice/cern.ch/user/z/zampolli/MergeCalib/OCDB

#ALIEN setting
# $1 = directory where to perform the find 
# $2 = run number
# $3 = OCDB path

echo ">>>>>>> Running AliRoot to merge calib objects found in $1"
aliroot -l -b -q merge.C\(\"$1\"\) 2>&1 | tee merge.log

echo ">>>>>>> Extract OCDB entries for run = $2, to be stored in $3"
aliroot -l -b -q makeOCDB.C\(\"$2\",\"$3\"\) 2>&1 | tee ocdb.log



