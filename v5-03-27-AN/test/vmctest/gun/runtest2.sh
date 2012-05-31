#!/bin/sh
# $Id$

# Script for reprocessing a given event number starting from saved random
# engine status files.
# Before running this script, you should run rungen.sh first.

NEVENTS=1
G4CONFIG="$ALICE_ROOT/test/vmctest/gun/g4Config2.C" 
RNDM_ROOT=$ALICE_ROOT/test/vmctest/gun/g4_ref/randomEvt4.root
RNDM_G4=$ALICE_ROOT/test/vmctest/gun/g4_ref/run0evt4.rndm
G4OUTDIR=g4

RUNG4=1

if [ "$RUNG4" = "1" ]; then 
  rm -rf *.root *.dat *.log fort* hlt hough raw* recraw/*.root recraw/*.log
  cp $RNDM_ROOT $RNDM_G4 .
  aliroot -b -q  sim.C\($NEVENTS,\""$G4CONFIG"\"\)  2>&1 | tee sim.log
  #aliroot -b -q rec.C      2>&1 | tee rec.log
  rm -fr $G4OUTDIR
  mkdir $G4OUTDIR
  mv *.root *.log *.rndm GRP $G4OUTDIR
  cp g4Config.C $G4OUTDIR
fi
