#!/bin/sh
# This script processes the ITS PID for n events
# (use the parameter n):
# - the AliTPCtest.C and AliITStest.C should be processed befor 
# (Dubnna version)
# 
# Do the PID procedure for the ITS. The output file AliITStrackV2Pid.root is
# created with the information for each track:
# the ITS ADC trancated mean signal, the particle reconstructed momentum,
# the probable weights for the different particle kinds (electron, pion, kaon,
# proton). These weights are under study now.
#
rm -f AliITStrackingV2Pid.root
#
aliroot -q -b "$ALICE_ROOT/ITS/save_pidV2.C(0,$[$1-1])" 
#
# This last line checks the ITS PID only, i.e. creates and draws
#the dEdx-momentum plot (comment if it's not necessary)
aliroot "$ALICE_ROOT/ITS/scan_pidV2.C(0,$[$1-1])"



