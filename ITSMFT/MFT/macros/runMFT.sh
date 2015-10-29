#!/bin/bash
date
# first declare default values
export NTRACKS=$1   # Number of tracks
NEVENTS=$2   # Number of events to be simulated
SIMCONFIG="Config.C"   # default simulation configuration file
RUN=169099  # run number for OCDB access
SEED=0  # random number generator seed should be used
SIMDIR="generated" # sub-directory where to move simulated files prior to reco
###############################################################################
# 
# Performing SIMULATION
# 
###############################################################################
  echo "Running simulation  ..."
  aliroot -l -b -q runSimulation.C\($NEVENTS,$RUN\) >testSim.out 2>testSim.err
  mkdir $SIMDIR
  echo "Copying generated files to $SIMDIR"
  cp Kinematics*.root galice.root TrackRefs*.root $SIMDIR
###############################################################################
# 
# Performing RECONSTRUCTION
# 
###############################################################################
  rm -f AliESD*.root *QA*.root
  echo "Running reconstruction  ..."
  aliroot -l -b -q runReconstruction.C >testReco.out 2>testReco.err

echo "Finished"  
ls -latr
