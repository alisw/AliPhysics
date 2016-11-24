#!/bin/bash -l
# The settings come from ~/.bash_profile

# Before running this script, you should run rungen.sh first.
# check if this is the case
if [ ! -f gen/galice.root ]; then
  echo "calling rungen.sh"
  ./rungen.sh
fi


NEVENTS=3
G3CONFIG="$ALICE_ROOT/test/vmctest/ppbench/g3Config.C" 
G4CONFIG="$ALICE_ROOT/test/vmctest/ppbench/g4Config.C" 
G3OUTDIR=g3
G4OUTDIR=g4

RUNG3=1
RUNG4=1

if [ "$RUNG3" = "1" ]; then 
  rm -rf *.root *.dat *.log fort* hlt hough raw* recraw/*.root recraw/*.log
  aliroot -b -q  sim.C\($NEVENTS,\""$G3CONFIG"\"\)  2>&1 | tee sim.log
  mv syswatch.log simwatch.log
  aliroot -b -q rec.C      2>&1 | tee rec.log
  mv syswatch.log recwatch.log
  aliroot -b -q ${ALICE_ROOT}/STEER/macros/CheckESD.C 2>&1 | tee check.log
  aliroot -b -q ${ALICE_ROOT}/STEER/macros/CreateAODfromESD.C 2>&1 | tee aod.log
  rm -fr $G3OUTDIR
  mkdir $G3OUTDIR
  mv *.root *.log *.ps GRP $G3OUTDIR
  cp g3Config.C $G3OUTDIR
fi

if [ "$RUNG4" = "1" ]; then 
  rm -rf *.root *.dat *.log fort* hlt hough raw* recraw/*.root recraw/*.log
  aliroot -b -q  sim.C\($NEVENTS,\""$G4CONFIG"\"\)  2>&1 | tee sim.log
  mv syswatch.log simwatch.log
  aliroot -b -q rec.C      2>&1 | tee rec.log
  mv syswatch.log recwatch.log
  aliroot -b -q ${ALICE_ROOT}/STEER/macros/CheckESD.C 2>&1 | tee check.log
  aliroot -b -q ${ALICE_ROOT}/STEER/macros/CreateAODfromESD.C 2>&1 | tee aod.log
  rm -fr $G4OUTDIR
  mkdir $G4OUTDIR
  mv *.root *.log *.rndm *.ps GRP $G4OUTDIR
  cp g4Config.C $G4OUTDIR
fi
