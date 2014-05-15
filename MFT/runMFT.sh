#!/bin/bash
date
# first declare default values
NEVENTS=10   # Number of events to be simulated
SIMCONFIG="Config.C"   # default simulation configuration file
CURDIR=`pwd`
echo 'current directory:' $CURDIR
OUTDIR=${CURDIR}/test
echo 'working directory:' $OUTDIR
mkdir $OUTDIR
RUN=169099  # run number for OCDB access
SIMDIR="generated" # sub-directory where to move simulated files prior to reco
# Copy *ALL* the macros we need in the output directory, not to mess
# with our source dir in any way.
cp $ALICE_ROOT/MFT/.rootrc \
   $ALICE_ROOT/MFT/Config.C \
   $ALICE_ROOT/MFT/rootlogon.C \
   $ALICE_ROOT/MFT/runReconstruction.C \
   $ALICE_ROOT/MFT/runSimulation.C \
   $ALICE_ROOT/MFT/AliMFTClusterQA.C \
   $ALICE_ROOT/MFT/AliMFTGeometry.root \
   $ALICE_ROOT/MFT/AODtrain.C \
   $ALICE_ROOT/MFT/RunAnalysisTaskMFTExample.C \
   $ALICE_ROOT/MFT/AliAnalysisTaskMFTExample.h \
   $ALICE_ROOT/MFT/AliAnalysisTaskMFTExample.cxx \
   $OUTDIR 
cd $OUTDIR

###############################################################################
# 
# Performing SIMULATION
# 
###############################################################################

echo "Running simulation  ..."
aliroot -l -b -q runSimulation.C\($NEVENTS,\""$SIMCONFIG"\",$RUN\) >$OUTDIR/testSim.out 2>$OUTDIR/testSim.err
mkdir $SIMDIR
echo "Copying generated files to $SIMDIR"
cp $OUTDIR/Kinematics*.root $OUTDIR/galice.root $OUTDIR/TrackRefs*.root $OUTDIR/$SIMDIR

###############################################################################
# 
# Performing RECONSTRUCTION and QA
# 
###############################################################################

rm -f AliESD*.root *QA*.root
echo "Running reconstruction  ..."
cd $OUTDIR

aliroot -l -b -q runReconstruction.C\(\""SAVEDIGITS"\"\) >$OUTDIR/testReco.out 2>$OUTDIR/testReco.err
      
aliroot -l -b -q AliMFTClusterQA.C\(\) >$OUTDIR/mftClusterQA.out 2>$OUTDIR/mftClusterQA.err

###############################################################################
# 
# Creating AODs
#
###############################################################################
  
echo "Creating AODs  ..."
echo aliroot -l -b -q AODtrain.C
aliroot -l -b -q AODtrain.C >$OUTDIR/AODLog.out 2>$OUTDIR/AODLog.err

###############################################################################
# 
# Performing analysis
#
###############################################################################
  
echo "Performing analysis  ..."
echo aliroot -l -b -q RunAnalysisTaskMFTExample.C
aliroot -l -b -q RunAnalysisTaskMFTExample.C >$OUTDIR/analysisLog.out 2>$OUTDIR/analysisLog.err

echo "Finished"  
echo "... see results in $OUTDIR"
ls -latr
tar czf out.tar.gz *.out
rm -f *.out *.so *.d
cd $CURDIR
