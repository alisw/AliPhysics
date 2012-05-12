#!/bin/bash

# Script to merge objects coming out of the calibration train:
# Arguments:
#    1 - directory on which to look for the files to be merged 
#    2 - runNumber
#    3 - OCDB output path
#    [4 - default OCDB]

# example:
# mergeMakeOCDB.sh /alice/cern.ch/user/a/aliprod/CPass0/output/ 120691 alien://folder=/alice/cern.ch/user/a/aliprod/CPass0/output

#ALIEN setting
# $1 = directory where to perform the find 
# $2 = runNumber
# $3 = OCDB path

path=$1
runNumber=$2
outputOCDB=$3

# if fourth argument given, its the default OCDB, otherwise use the default raw://
defaultOCDB="raw://"
[[ -n $4 ]] && defaultOCDB=$4

# if fifth argument given it is a runnuber, otherwise guess from path (when ran on alien)
runNumber=$(echo $1 | cut -d "/" -f 6 | sed 's/^0*//')
[[ -n $5 ]] && runNumber=$5

if [ -f Run0_999999999_v3_s0.root ]; then
    mkdir -p TPC/Calib/Correction
    mv Run0_999999999_v3_s0.root TPC/Calib/Correction/
fi

echo ">>>>>>> Running AliRoot to merge calib objects found in $path with pattern AliESDfriends_v1.root"
aliroot -l -b -q merge.C\(\"$path\",\"AliESDfriends_v1.root\"\) 2>&1 | tee merge.log
mv syswatch.log syswatch_merge.log

echo ">>>>>>> Extract OCDB entries for run = $runNumber, to be stored in $outputOCDB"
aliroot -l -b -q makeOCDB.C\($runNumber,\"$outputOCDB\",\"$defaultOCDB\"\) 2>&1 | tee ocdb.log
mv syswatch.log syswatch_makeOCDB.log
