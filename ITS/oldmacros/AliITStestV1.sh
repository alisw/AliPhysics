#!/bin/sh

# delete eventual old files from the last run
./AliITSDeleteOldFiles.sh

# do everything for the TPC
aliroot -q -b "$ALICE_ROOT/TPC/AliTPCtest.C"

# create summable digits for the ITS
aliroot -q -b "$ALICE_ROOT/ITS/AliITSHits2SDigits.C"

# go from summable digits to digits for the ITS
aliroot -q -b "$ALICE_ROOT/ITS/AliITSSDigits2Digits.C"

# create reconstructed points for the ITS
aliroot -q -b "$ALICE_ROOT/ITS/AliITSDigits2RecPoints.C"

# do the tracking V1
aliroot -q -b "$ALICE_ROOT/ITS/AliITSTrackingV1.C"

# prepare results of tracking V1 for comparison
aliroot -q -b "$ALICE_ROOT/ITS/AliITSStoreFindableTracks.C"

# do ITS tracking V1 comparison
aliroot -q -b "$ALICE_ROOT/ITS/AliITSComparisonV1.C"

#
# after all of the above you can run AliITSPlotTracksV1.C macro under 
# aliroot to see plots of the efficiency and track parameter resolution
# or AliITSDisplayTracksV1.C to display found tracks on top of the
# ITS detailed geometry
