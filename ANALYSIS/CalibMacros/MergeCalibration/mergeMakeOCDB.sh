#!/bin/bash

# Script to merge objects coming out of the calibration train:
# Arguments:
#    1 - directory on which to look for the files to be merged 
#    2 - OCDB output path
# example:
# mergeMakeOCDB.sh /alice/cern.ch/user/z/zampolli/CalibTrain/output/

#ALIEN setting
# $1 = directory where to perform the find 
#runnum=$1
#ocdbStorage=$2
#echo Directory to look into = $1

echo Run path $1
echo OCDB output path $2
runNum=`basename $1`
echo RunNumber $runNum

echo ">>>>>>> Running AliRoot to merge calib objects found in $1"
aliroot -l -b -q merge.C\(\"$1\"\) 2>&1 | tee merge.log

echo ">>>>>>> Extract OCDB entries"
aliroot -l -b -q makeOCDB.C\(\"$runNum\",\"$2\"\) 2>&1 | tee ocdb.log



