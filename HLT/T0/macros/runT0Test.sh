#!/bin/bash

# -------------------------------------------
# Test full HLT chains with TZERO reconstruction
# -------------------------------------------

# N events
NEVENTS=20

inputFile=${1}
firstEvent=1
lastEvent=5

aliroot -l -q -b testT0Config.C $ALICE_ROOT/HLT/exa/recraw-local.C"(\"$inputFile\",\"local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2016/OCDB/\", $firstEvent, $lastEvent, \"HLT\", \"chains=RootWriter ignore-hltout TPC-input=compressed -loglevel=0x7c\")"
