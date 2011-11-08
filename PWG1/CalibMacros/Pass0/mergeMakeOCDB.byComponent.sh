#! /bin/bash

# init
path=$1
run=$2
ocdb=$3
echo "***********************" 2>&1 | tee -a merge.log
echo mergeMakeOCDB.sh started 2>&1 | tee -a merge.log
echo path = $path 2>&1 | tee -a merge.log
echo run  = $run 2>&1 | tee -a merge.log
echo ocdb = $ocdb 2>&1 | tee -a merge.log
echo "***********************" 2>&1 | tee -a merge.log

# setup components
components="TOF MeanVertex T0 TRD TPC SDD"

# copy
echo "***********************" 2>&1 | tee -a merge.log
echo copying files for run $run 2>&1 | tee -a merge.log
echo from $path 2>&1 | tee -a merge.log
echo "***********************" 2>&1 | tee -a merge.log
aliroot -b -q "merge.C(\"COPY\", \"$path\")" 2>&1 | tee -a merge.log
mv syswatch.log copy_syswatch.log

# process by component
for det in $components; do

    # merge
    echo "***********************" 2>&1 | tee -a merge.log
    echo merging $det data 2>&1 | tee -a merge.log
    echo "***********************" 2>&1 | tee -a merge.log
    aliroot -b -q "merge.C(\"$det\", \"calib.list\")" 2>&1 | tee -a merge.log
    mv syswatch.log $det\_merge_syswatch.log

    # make OCDB
    echo "***********************" 2>&1 | tee -a ocdb.log
    echo making $det OCDB 2>&1 | tee -a ocdb.log
    echo "***********************" 2>&1 | tee -a ocdb.log
    aliroot -b -q "makeOCDB.C(\"CalibObjects.root\", \"$det\", \"$run\", \"$ocdb\")" 2>&1 | tee -a ocdb.log
    mv CalibObjects.root $det\_CalibObjects.root

done

# global merge
echo "***********************" 2>&1 | tee -a merge.log
echo merging ALL data 2>&1 | tee -a merge.log
echo "***********************" 2>&1 | tee -a merge.log
ls *CalibObjects.root > objects.list
aliroot -b -q "merge.C(\"ALL\", \"objects.list\")" 2>&1 | tee -a merge.log
touch CalibObjects.root

# summary
echo "***********************" 2>&1 | tee -a ocdb.log
echo SUMMARY 2>&1 | tee -a ocdb.log
echo "***********************" 2>&1 | tee -a ocdb.log
ls -altr *CalibObjects.root *done 2>&1 | tee -a ocdb.log