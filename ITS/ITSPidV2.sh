#!/bin/sh
# delete eventual old files from the last run
./DeleteOldV2Files
# run the hit generation
aliroot -q -b "$ALICE_ROOT/macros/grun.C"
# digitize TPC and prepare TPC tracks for matching with the ITS
aliroot -q -b "$ALICE_ROOT/TPC/AliTPCtest.C"
# digitize ITS
aliroot -q -b "$ALICE_ROOT/ITS/AliITSHitsToDigitsDefault.C"
# create reconstructed points for the ITS
aliroot -q -b "$ALICE_ROOT/ITS/ITSDigitsToClusters.C"
#prepare slow or fast recpoints for the tracking 
aliroot -q -b "$ALICE_ROOT/ITS/"AliITSFindClustersV2.C\(\'$2\'\)
# ITS tracking
aliroot -q -b "$ALICE_ROOT/ITS/SetConvConst.C"  "$ALICE_ROOT/ITS/AliITSFindTracksV2.C"
# do the PID
aliroot -q -b "$ALICE_ROOT/ITS/SetConvConst.C" "$ALICE_ROOT/ITS/save_pid.C"
aliroot "$ALICE_ROOT/ITS/ITSPIDtest.C"
#
# After all of the above the PID menu will appear. To check the PID package  
# put the buttons:
#
#1.  "fill tab_tr" and "dEdX-P plot" (or dEdX-P pions, ... kaons, ... elect, 
# ... prot) and the dEdX versus momentum taken from the tracking for all
# particle (pions, kaons, electrons, protons) will be drown.
#
#2. "dEdX.C" and the same other ones as for the point 1.  In this case the 
# same plots but fot the generated particle momentum will be drown. 
#
# The output file AliITStrackV2Pid.root will be written with the information
# for each track:
# the ITS ADC trancated mean signal, the particle reconstructed momentum,
# the probable weights for the different particle kinds (electron, pion, kaon,
# proton). These weights are under study now.








