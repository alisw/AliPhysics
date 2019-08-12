#! /bin/bash

inputDir=$1
lhcPeriod=$2
runListFile=$3
automaticMode=$4
nProcs=$5

macrosPath=$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/calibration/

mkdir -p temp/anchors
for run in `cat $runListFile`; do

   run=${run//[^0-9]}
   
   if [ -f temp/anchors/Anchor_${lhcPeriod}_${run}_VHM.txt ]; then
      echo "Anchor already existing for run $run"
      continue
   fi

   if [ $automaticMode == 'yes' ]
   then
      while [ $(ps -ef | grep aliroot | grep DetermineAnchorsPP | wc -l) -ge $nProcs ]; do
         sleep 5s;
      done
      
      nohup aliroot -b -q "$macrosPath/DetermineAnchorsPP.C(\"$inputDir\", \"$lhcPeriod\",$run,\"\", kTRUE)" &> temp/nohupOutput/nohupDetermineAnchors_$run.txt &
      
   else
      aliroot -q "$macrosPath/DetermineAnchorsPP.C(\"$inputDir\", \"$lhcPeriod\",$run,\"\", kFALSE)"
   fi

done
