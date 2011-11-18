#!/bin/sh
# $Id $

NEVENTS=1
G3CONFIG="$ALICE_ROOT/test/vmctest/lego/g3Config.C" 
G4CONFIG="$ALICE_ROOT/test/vmctest/lego/g4Config.C" 

RUNG3=1
RUNG4=1

#for DET in ABSO DIPO FMD FRAME HALL ITS MAG MUON PHOS PIPE PMD HMPID SHIL T0 TOF TPC TRD ZDC EMCAL ACORDE VZERO; do
for DET in ALL; do

if [ "$RUNG3" = "1" ]; then 
  G3OUTDIR=g3/lego_$DET
  rm -rf *.root *.dat *.log fort* hlt hough raw* recraw/*.root recraw/*.log
  aliroot -b -q  sim.C\($NEVENTS,\"$G3CONFIG\",\"$DET\"\)  2>&1 | tee sim.log
  rm -fr $G3OUTDIR
  mkdir -p $G3OUTDIR
  mv *.root *.log $G3OUTDIR
  cp $G3CONFIG $G3OUTDIR
fi

if [ "$RUNG4" = "1" ]; then 
  G4OUTDIR=g4/lego_$DET
  rm -rf *.root *.dat *.log fort* hlt hough raw* recraw/*.root recraw/*.log
  aliroot -b -q  sim.C\($NEVENTS,\"$G4CONFIG\",\"$DET\"\)  2>&1 | tee sim.log
  rm -fr $G4OUTDIR
  mkdir -p $G4OUTDIR
  mv *.root *.log $G4OUTDIR
  cp $G4CONFIG $G4OUTDIR
fi

done
