#!/bin/sh
# This script processes the following steps for n events
# (use the parameter n):
# - the TPC simulation and tracking (Marian Ivanov version),
# - the TPC PID,
# - the AliTPCtracksPid.root is created with the following information:
#    fLabel - track number
#    fPcode - particle code after the TPC PID
#    fMom   - particle momentum      (from the AliTPCtracks.root)
#    fLam   - particle lambder angle (from the AliTPCtracks.root)
#    fPhi   - particle phi angle     (from the AliTPCtracks.root)    
#    fSignal- TPC signal (mips)
#    fWpi   - the PID weight for the identified pion
#    fWk    - the PID weight for the identified kaon
#    fWp    - the PID weight for the identified proton
# and the AliTPCScanPID.C is the example of macros to read the PID
# information. 
# See PIDplot.ps after run of this script.
# (Dubna , Aug 2002, Jan 2003)
# 
if [ $# = 0   ]; then nev=1; else nev=$1; fi
if [ $nev = 0 ]; then nev=1; fi
# delete eventual old files from the last run
echo "Start simulation for " $nev " event(s)"
rm -f *.root PIDplot.ps
#
# run the hit generation
aliroot -q -b "$ALICE_ROOT/macros/grun.C($nev)" 
#
# TPC simulation 
aliroot -q -b "$ALICE_ROOT/TPC/AliTPCHits2Digits.C($nev)"

# TPC tracking (MI version)
aliroot -q -b "$ALICE_ROOT/TPC/AliTPCFindClustersMI.C($nev)" 
aliroot -q -b "$ALICE_ROOT/TPC/AliTPCFindTracksMI.C($nev)" 
#
# PID in TPC
aliroot -q -b "AliTPCSavePID.C($nev)"
aliroot -q -b "AliTPCScanPID.C($nev)"
#


