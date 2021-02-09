#! /bin/bash

inputDirectory=$1
lhcPeriod=$2
defaultRun=$3
nProcs=$4
runListFile=$5

macrosPath=$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/

mkdir -p temp/nohupOutput
mkdir -p temp/buffers
mkdir -p temp/partialOADBs
files=`ls $inputDirectory`

echo "Running the calibration for MB triggers..."
for ifile in $files; do
   # The format of the filenames should be: AnalysisResults_<identifier>.root, where "<identifier>" is either just the run number or some other name in case
   # this file is obtained from merging more runs.

   # find out the file identifier (run number or whatever)
   runIdentifier=${ifile#AnalysisResults_}
   runIdentifier=${runIdentifier%.root}
    
   # limit the number of simultaneous processes
   while [ $(ps -ef | grep aliroot | grep CalibratePeriodPP | wc -l) -ge $nProcs ]; do
      sleep 5s;
   done
    
   nohup aliroot –l –b –q "$macrosPath/calibration/CalibratePeriodPP.C(\"$inputDirectory\", \"$lhcPeriod\", \"MB\", $defaultRun, \"$runIdentifier\")" &> temp/nohupOutput/nohupCalibMB_$runIdentifier.txt &
done

# wait for all the calibration jobs to finish
while [ $(ps -ef | grep aliroot | grep CalibratePeriodPP | wc -l) -ge 1 ]; do
   sleep 10s;
done

# clean the buffer files
rm -r temp/buffers/*.root

# merge the individual OADB files 
echo "Merging individual OADB files..."
ls temp/partialOADBs/OADB-$lhcPeriod*MB.root > oadbFiles.txt
cat oadbFiles.txt
aliroot -l -b -q "$macrosPath/calibration/MergeOADB.C(\"oadbFiles.txt\",\"OADB-$lhcPeriod-MB.root\")"
rm oadbFiles.txt

# run the OADB test for every run and produce pictures
echo "Run the OADB tests..."
mkdir -p temp/oadbTests
for run in `cat $runListFile`; do
 
   run=${run//[^0-9]}
 
   while [ $(ps -ef | grep aliroot | grep TestOADBMultSelPP | wc -l) -ge $nProcs ]; do
      sleep 5s;
   done
   
   nohup aliroot –b –q "$macrosPath/qa/TestOADBMultSelPP.C(\"$inputDirectory\", \"$lhcPeriod\", $run)" &> temp/nohupOutput/nohupTestOADB_$run.txt &
done

while [ $(ps -ef | grep aliroot | grep TestOADBMultSelPP | wc -l) -ge 1 ]; do
   sleep 10s;
done

echo "Done!"
