#!/bin/sh
# 9 July 2002,Dubna
# This script processes the following steps for n events (AliRoot v3-08-03)
# (use the parameter n):
# - the TPC and ITS slow simulation (hits+digits+clusters),
# - the TPC+ITS tracking,
# - the TPC and ITS PID
# (Dubnna version)
# 
if [ $# = 0   ]; then nev=1; else nev=$1; fi
if [ $nev = 0 ]; then nev=1; fi
# delete eventual old files from the last run
echo "Start simulation for " $nev " event(s)"
$ALICE_ROOT/ITS/AliITSDeleteOldFiles.sh
#
# run the hit generation
aliroot -q -b "$ALICE_ROOT/macros/grun.C($nev)" 
#
# digitize TPC
aliroot -q -b "$ALICE_ROOT/TPC/AliTPCHits2Digits.C($nev)"
#
# digitize ITS
aliroot -q -b "$ALICE_ROOT/ITS/AliITSMakeSDigits.C($nev)"
aliroot -q -b "$ALICE_ROOT/ITS/AliITSSDigits2Digits.C"
# create reconstructed points for the ITS
aliroot -q -b "$ALICE_ROOT/ITS/AliITSDigits2RecPoints.C"

# do the TPC+ITS tracking

#bad: aliroot -q -b "$ALICE_ROOT/ITS/AliBarrelReconstructionV2.C($1)" 
aliroot -q -b "$ALICE_ROOT/ITS/AliBarrelReconstructionV2.C($nev)"
#
# Do the PID procedure for the ITS. The output file AliITStrackV2Pid.root is
# created with the information for each track:
# the ITS ADC trancated mean signal, the particle reconstructed momentum,
# the probable weights for the different particle kinds (electron, pion, kaon,
# proton). These weights are under study now.
#
aliroot -q -b "$ALICE_ROOT/ITS/AliITSSavePIDV2.C(0,$[$nev-1])" 
#
# This last line checks the ITS PID only, i.e. creates and draws
#the dEdx-momentum plot (comment if it's not necessary)
aliroot "$ALICE_ROOT/ITS/AliITSScanPIDV2.C(0,$[$nev-1])"


