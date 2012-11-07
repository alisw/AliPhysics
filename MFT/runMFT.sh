#!/bin/bash
date
# first declare default values
NEVENTS=10   # Number of events to be simulated
SIMCONFIG="Config.C"   # default simulation configuration file
CURDIR=`pwd`
echo $CURDIR
OUTDIR=${CURDIR}/test
echo $OUTDIR
mkdir $OUTDIR
RUN=169099  # run number for OCDB access
SEED=12345  # random number generator seed should be used
SIMDIR="generated" # sub-directory where to move simulated files prior to reco
# Copy *ALL* the macros we need in the output directory, not to mess
# with our source dir in any way.
cp $ALICE_ROOT/MFT/.rootrc \
   $ALICE_ROOT/MFT/Config.C \
   $ALICE_ROOT/MFT/rootlogon.C \
   $ALICE_ROOT/MFT/runReconstruction.C \
   $ALICE_ROOT/MFT/runSimulation.C \
   $ALICE_ROOT/MFT/AliMuonForwardTrackFinder.C \
   $ALICE_ROOT/MFT/FilterMuonGlobalTracks.C \
   $ALICE_ROOT/MFT/AliMFTClusterQA.C \
   $ALICE_ROOT/MFT/AliMFTGeometry.root \
   $OUTDIR 
cd $OUTDIR
###############################################################################
# 
# Performing SIMULATION
# 
###############################################################################
  echo "Running simulation  ..."
  aliroot -l -b -q runSimulation.C\($SEED,$NEVENTS,\""$SIMCONFIG"\",$RUN\) >$OUTDIR/testSim.out 2>$OUTDIR/testSim.err
  mkdir $SIMDIR
  echo "Copying generated files to $SIMDIR"
  cp $OUTDIR/Kinematics*.root $OUTDIR/galice.root $OUTDIR/TrackRefs*.root $OUTDIR/$SIMDIR
###############################################################################
# 
# Performing RECONSTRUCTION
# 
###############################################################################
  rm -f AliESD*.root *QA*.root
  echo "Running reconstruction  ..."
  cd $OUTDIR

  aliroot -l -b -q runReconstruction.C\($SEED,\""SAVEDIGITS"\"\) >$OUTDIR/testReco.out 2>$OUTDIR/testReco.err

  aliroot -l -b -q AliMFTClusterQA.C\(\) >$OUTDIR/mftClusterQA.out 2>$OUTDIR/mftClusterQA.err

  aliroot -l -b -q AliMuonForwardTrackFinder.C\($RUN,1\) >$OUTDIR/globalTracking.withBransonCorrection.out 2>$OUTDIR/globalTracking.withBransonCorrection.err
  mv MuonGlobalTracking.QA.run$RUN.root MuonGlobalTracking.QA.run$RUN.withBransonCorrection.root 
  mv MuonGlobalTracks.root MuonGlobalTracks.withBransonCorrection.root 

  aliroot -l -b -q AliMuonForwardTrackFinder.C\($RUN,0\) >$OUTDIR/globalTracking.withoutBransonCorrection.out 2>$OUTDIR/globalTracking.withoutBransonCorrection.err
  mv MuonGlobalTracking.QA.run$RUN.root MuonGlobalTracking.QA.run$RUN.withoutBransonCorrection.root 
  mv MuonGlobalTracks.root MuonGlobalTracks.withoutBransonCorrection.root 

  aliroot -l -b -q FilterMuonGlobalTracks.C+\(\) >$OUTDIR/filterGlobalTracks.out 2>$OUTDIR/filterGlobalTracks.err

echo "Finished"  
echo "... see results in $OUTDIR"
ls -latr
tar czf out.tar.gz *.out
rm -f *.out *.so *.d
cd $CURDIR
