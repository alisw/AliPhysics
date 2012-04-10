#!/bin/bash

# Script to merge objects coming out of the calibration train:
# Arguments:
#    1 - directory on which to look for the files to be merged 
#    2 - pattern
#    3 - OCDB output path

# example:
# mergeMakeOCDB.sh /alice/cern.ch/user/a/aliprod/CPass1/output/ 000120691 alien://folder=/alice/cern.ch/user/a/aliprod/CPass1/output

#ALIEN setting
# $1 = directory where to perform the find 
# $2 = pattern
# $3 = OCDB path

if [ -f Run0_999999999_v3_s0.root ]; then
    mkdir -p TPC/Calib/Correction
    mv Run0_999999999_v3_s0.root TPC/Calib/Correction/
fi

echo ">>>>>>> Running AliRoot to merge calib objects found in $1 with pattern AliESDfriends_v1.root"
aliroot -l -b -q merge.C\(\"$1\",\"AliESDfriends_v1.root\"\) 2>&1 | tee merge.log
mv syswatch.log syswatch_merge.log

echo ">>>>>>> Extract OCDB entries for run = $2, to be stored in $3"
aliroot -l -b -q makeOCDB.C\(\"$2\",\"$3\"\) 2>&1 | tee ocdb.log
mv syswatch.log syswatch_makeOCDB.log
