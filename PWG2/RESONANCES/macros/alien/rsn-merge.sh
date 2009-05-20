#!/bin/sh
#
# Calls the merging macro and passes to it the two required arguments
# -- $1 = root path containing all sub-paths with files to merge
# -- $2 = name of the files to be merged
# -- $3 = name of output file

path=/alice/cern.ch/user/p/pulvir/$1
outName=$3.root
logName=$3.log
inName=$2.root

exec root -b -q -l -x RsnMergeAlien.C --path "$path" --name "$inName" --out "$outName" >& $logName
echo >> $logName
echo >> $logName
echo -n 'Executed in: ' >> $logName
hostname -f >> $logName
