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
runNumber=`echo "$2" | sed 's/^0*//'`
outputOCDB=$3

# if fourth argument given, it is the default OCDB, otherwise use the default raw://
defaultOCDB="raw://"
[[ -n $4 ]] && defaultOCDB=$4

if [ -f Run0_999999999_v3_s0.root ]; then
    mkdir -p TPC/Calib/Correction
    mv Run0_999999999_v3_s0.root TPC/Calib/Correction/
fi

echo "* Running AliRoot to merge calib objects found in $path with pattern CalibObjects.root"
echo aliroot -l -b -q "merge.C(\"$path\",\"CalibObjects.root\")"
time aliroot -l -b -q "merge.C(\"$path\",\"CalibObjects.root\")" &> merge.log

exitcode=$?

echo "*! Exit code: $exitcode"

mv syswatch.log syswatch_merge.log

echo "* Extract OCDB entries for run = $runNumber, to be stored in $outputOCDB"
echo aliroot -l -b -q "makeOCDB.C($runNumber,\"$outputOCDB\",\"$defaultOCDB\")"
time aliroot -l -b -q "makeOCDB.C($runNumber,\"$outputOCDB\",\"$defaultOCDB\")" &> ocdb.log

exitcode=$?

echo "*! Exit code: $exitcode"

mv syswatch.log syswatch_makeOCDB.log
