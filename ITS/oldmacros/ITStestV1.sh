#!/bin/sh
# delete eventual old files from the last run
./DeleteOldFiles
# run the hit generation
aliroot -q -b "$ALICE_ROOT/macros/grun.C($1)"
# digitize TPC
aliroot -q -b "$ALICE_ROOT/TPC/AliTPCHits2Digits.C($1)"
# prepare TPC tracks for matching with the ITS
aliroot -q -b "$ALICE_ROOT/ITS/TPCtracking.C($1)"
# digitize ITS
aliroot -q -b "$ALICE_ROOT/ITS/AliITSHits2DigitsDefault.C"
# create reconstructed points for the ITS
aliroot -q -b "$ALICE_ROOT/ITS/ITSDigitsToClusters.C(0,$[$1-1])"
# prepare for tracking
aliroot -q -b "$ALICE_ROOT/ITS/ITStracks.C(0,$[$1-1])"
# do the tracking
aliroot -q -b "$ALICE_ROOT/ITS/ITStracking.C(0,$[$1-1])"
# do the comparison
aliroot -q -b "$ALICE_ROOT/ITS/AliITSComparisonV1.C(0,$[$1-1])"
#
# after all of the above you can run ITSPlotTracks.C macro under aliroot to
# see plots of the efficiency and track parameter resolution
