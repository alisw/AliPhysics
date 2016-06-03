#!/bin/bash

dirname=
pattern="ResidualTrees.root"

[[ $# < 2 ]] &&  echo 'procRunRes.sh <path> <run>' &&  exit

dirname=$1
run=$2
inpList="lst.txt"

echo "...searching for data as:"
echo alien_find "$dirname" "$pattern" '| grep ' "$dirname" ' | perl -p -e ' 's/$dirname/alien:\/\/$dirname/'
alien_find "$dirname" "$pattern" | grep "$dirname" | perl -p -e 's/$dirname/alien:\/\/$dirname/' >& $inpList
chunks=`wc -l $inpList | cut -f1 -d' '`

[[ $chunks < 1 ]] && echo "did not find any $pattern for run $run in $dirname" && exit

curdir=`pwd`

scriptVDT="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/procVDTime.sh"
scriptMAP="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/procRunTimeBins.sh"

[[ ! -f "$scriptVDT" ]] && echo "did not find script $scriptVDT" && exit
[[ ! -f "$scriptMAP" ]] && echo "did not find script $scriptMAP" && exit 

runDir="$curdir""/resmap""$run"
echo "will process run $run in $runDir"
[ ! -d $runDir  ] && mkdir $runDir
mv $inpList $runDir
cd $runDir

echo "Running vdrift extraction"
$scriptVDT "$inpList" "$run"
echo "Running map extraction"
$scriptMAP "timeBins.log"

cd $curDir

exit 0

