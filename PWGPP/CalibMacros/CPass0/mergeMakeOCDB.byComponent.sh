#! /bin/bash

# Script to merge objects coming out of the calibration train and export of the OCDB:
# Arguments:
#    1 - directory on which to look for the files to be merged 
#    2 - pattern
#    3 - OCDB output path
#    [4 - input default OCDB]
#

#ALIEN setting
# $1 = directory where to perform the find 
# $2 = pattern
# $3 = OCDB path

# if fourth argument given, its the default OCDB, otherwise use the default raw://
defaultOCDB="raw://"
[[ $# -eq 4 ]] && defaultOCDB=$4

# init
path=$1
run=$2
ocdb=$3
isLocal=0
[[ -f $path ]] && isLocal=1

echo "***********************" 2>&1 | tee -a merge.log
echo mergeMakeOCDB.byComponent.sh started 2>&1 | tee -a merge.log
echo path = $path 2>&1 | tee -a merge.log
echo run  = $run 2>&1 | tee -a merge.log
echo ocdb = $ocdb 2>&1 | tee -a merge.log
echo isLocal = $isLocal 2>&1 | tee -a merge.log
echo "***********************" 2>&1 | tee -a merge.log

# setup components
components="TOF MeanVertex T0 SDD TRD TPCCalib  TPCAlign TPCCluster"

# copy
if [ $isLocal -eq 0 ]; then
    echo "***********************" 2>&1 | tee -a merge.log
    echo copying files for run $run 2>&1 | tee -a merge.log
    echo from $path 2>&1 | tee -a merge.log
    echo "***********************" 2>&1 | tee -a merge.log
    aliroot -b -q "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/mergeByComponent.C(\"COPY\",0,  \"$path\", \"AliESDfriends_v1.root\")" 2>&1 | tee -a merge.log
#mv syswatch.log copy_syswatch.log
else
  cp $path calib.list
fi;

# process by component
for det in $components; do

    # merge
    echo "***********************" 2>&1 | tee -a merge.log
    echo merging $det data 2>&1 | tee -a merge.log
    echo "***********************" 2>&1 | tee -a merge.log
    aliroot -b -q "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/mergeByComponent.C(\"$det\", \"calib.list\")" 2>&1 | tee -a merge.log
    mv syswatch.log merge_syswatch_$det.log
    mv CalibObjects.root CalibObjects_$det.root
done


# global merge
echo "***********************" 2>&1 | tee -a merge.log
echo merging ALL data 2>&1 | tee -a merge.log
echo "***********************" 2>&1 | tee -a merge.log
ls CalibObjects*.root > objects.list
aliroot -b -q "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/mergeByComponent.C(\"ALL\", \"objects.list\")" 2>&1 | tee -a merge.log
#touch CalibObjects.root


# make OCDB
echo "***********************" 2>&1 | tee -a ocdb.log
echo making $det OCDB 2>&1 | tee -a ocdb.log
echo "***********************" 2>&1 | tee -a ocdb.log
aliroot -b -q "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/makeOCDB.C(\"$run\", \"$ocdb\", \"$defaultOCDB\")" 2>&1 | tee -a ocdb.log
#mv CalibObjects.root $det\_CalibObjects.root


# summary
echo "***********************" 2>&1 | tee -a ocdb.log
echo SUMMARY 2>&1 | tee -a ocdb.log
echo "***********************" 2>&1 | tee -a ocdb.log
ls -altr *CalibObjects.root *done 2>&1 | tee -a ocdb.log