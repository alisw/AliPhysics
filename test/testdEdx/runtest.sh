#!/bin/sh

#
# Here we should write waht is the aim of test
#

export TestdEdxNTracks=$1
nevents=$2
echo Running dEdx digitzer test job 
echo NEvents $nevents
echo NTracks per event  $TestdEdxNTracks


rm -rf *.root *.dat *.log fort* hlt hough raw* recraw/*.root recraw/*.log
echo aliroot -b -q sim.C\($nevents\)     
aliroot -b -q sim.C\($nevents\)      2>&1 | tee sim.log
mv syswatch.log simwatch.log
aliroot -b -q $1rec.C      2>&1 | tee rec.log
mv syswatch.log recwatch.log
aliroot -b -q ${ALICE_ROOT}/STEER/CheckESD.C 2>&1 | tee check.log
aliroot -b -q aod.C 2>&1 | tee aod.log

cd recraw
ln -s ../raw.root .
aliroot -b -q rec.C      2>&1 | tee rec.log
mv syswatch.log ../rawwatch.log
aliroot -b -q aod.C 2>&1 | tee aod.log

exit;


submitMultiplicityScan(){
#
# Here we submit the jobs for the simulation//reconstruction
# Parameters:
#   nbins to be investigated   (default 5)
#   max multiplicity per event (default central PbPb - 15000)
#   minimal total statistic of the tracks (default 4 central events - 60000 tracks)
# Jobs will be submitted per event     
#
# For each setting new directory will be created - indicating muiltiplicity
# dir<ntracks>/dir<eventNr>  

  
}
