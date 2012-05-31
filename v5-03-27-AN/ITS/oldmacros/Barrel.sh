#!/bin/sh
# This script processes the following steps for n events
# (use the parameter n):
# - the TPC and ITS slow simulation (hits+digits+clusters),
# - the TPC+ITS tracking,
# - the TPC and ITS PID
# (Dubnna version)
# 
# delete eventual old files from the last run
rm -f *.root
#
# run the hit generation
aliroot -q -b "$ALICE_ROOT/macros/grun.C($1)" 
#
# digitize TPC
aliroot -q -b "$ALICE_ROOT/TPC/AliTPCHits2Digits.C($1)"
#
# digitize ITS
aliroot -q -b "$ALICE_ROOT/ITS/AliITSHits2DigitsDefault.C" 
#
# create reconstructed points for the ITS
aliroot -q -b "$ALICE_ROOT/ITS/ITSDigitsToClusters.C(0,$[$1-1])" 
#
# do the TPC+ITS tracking
aliroot -q -b "$ALICE_ROOT/ITS/AliBarrelReconstructionV2.C($1)" 
#
# Do the PID procedure for the ITS. The output file AliITStrackV2Pid.root is
# created with the information for each track:
# the ITS ADC trancated mean signal, the particle reconstructed momentum,
# the probable weights for the different particle kinds (electron, pion, kaon,
# proton). These weights are under study now.
#
aliroot -q -b "$ALICE_ROOT/ITS/save_pidV2.C(0,$[$1-1])" 
#
# This last line checks the ITS PID only, i.e. creates and draws
#the dEdx-momentum plot (comment if it's not necessary)
aliroot "$ALICE_ROOT/ITS/scan_pidV2.C(0,$[$1-1])"


