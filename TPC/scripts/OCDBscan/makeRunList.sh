#!/bin/bash

# Create the file  run.list to be processed
# Aruments
# 1   -  prefix
# Example usage
# $ALICE_ROOT/TPC/scripts/OCDBscan/makeRunList.sh /alice/data/2010


prefix=$1
alien_find $prefix/OCDB/GRP/GRP/Data Run            > grp.list
alien_find $prefix/OCDB/TPC/Calib/HighVoltage  Run  > hv.list 

for afile in `cat hv.list | grep root`; do 
    bfile=`basename $afile`; 
    echo $bfile  | sed s/_/\ / | sed s_Run__ | gawk '{print $1 }' 
done &> hvRun.list
sort  hvRun.list > hv.list

for arun  in `cat hv.list`; do 
    grpstatus=`cat grp.list | grep $arun`
    if [ -n $grpstatus ] ; then
      echo $arun
    fi;
done > hvRun.list

sort  hvRun.list > run.list

