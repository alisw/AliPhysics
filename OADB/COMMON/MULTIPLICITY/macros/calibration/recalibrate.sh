#! /bin/bash

inputDir=$1
lhcPeriod=$2
defaultRun=$3
nProcs=$4
runListFile=$5


macrosPath=$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/

mkdir -p temp/partialOADBs
for run in `cat $runListFile`; do

   run=${run//[^0-9]}
      
   while [ $(ps -ef | grep aliroot | grep ReCalibratePeriodPP | wc -l) -ge $nProcs ]; do
      sleep 5s;
   done
      
   nohup aliroot -b -q "$macrosPath/calibration/ReCalibratePeriodPP.C(\"$inputDir\", $run, \"$lhcPeriod\", \"VHM\", $defaultRun)" &> temp/nohupOutput/nohupReCalib_$run.txt &

done

# wait for all the calibration jobs to finish
while [ $(ps -ef | grep aliroot | grep ReCalibratePeriodPP | wc -l) -ge 1 ]; do
   sleep 10s;
done

# clean the buffer files
rm -r temp/buffers/*.root

# merge the individual OADB files 
echo "Merging individual OADB files..."
ls temp/partialOADBs/OADB-$lhcPeriod*VHM.root > oadbFiles.txt
cat oadbFiles.txt
aliroot -l -b -q "$macrosPath/calibration/MergeOADB.C(\"oadbFiles.txt\",\"OADB-$lhcPeriod-VHM.root\")"
rm oadbFiles.txt

# stich the MB and VHM OADBs 
aliroot -l -b -q "$macrosPath/calibration/StitchOADBs.C(\"$lhcPeriod\", $defaultRun)"

# test stitching
mkdir -p temp/checkStitching
for run in `cat $runListFile`; do

   run=${run//[^0-9]}
      
   while [ $(ps -ef | grep aliroot | grep TestStitchedOADB | wc -l) -ge $nProcs ]; do
      sleep 5s;
   done
      
   nohup aliroot -b -q "$macrosPath/qa/TestStitchedOADB.C(\"$inputDir\", \"$lhcPeriod\", $run)" &> temp/nohupOutput/nohupTestStitch_$run.txt &

done
