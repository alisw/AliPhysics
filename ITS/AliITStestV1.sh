#!/bin/sh

# delete eventual old files from the last run
./AliITSDeleteOldFiles.sh

# run the hit generation
aliroot -q -b "$ALICE_ROOT/macros/grun.C($1)"

# digitize TPC
aliroot -q -b "$ALICE_ROOT/TPC/AliTPCHits2Digits.C($1)"

# prepare TPC tracks for matching with the ITS
aliroot -q -b "$ALICE_ROOT/ITS/AliTPCTracking4ITS.C($1)"

# create summable digits for the ITS
aliroot -q -b "$ALICE_ROOT/ITS/AliITSHits2SDigits.C"

# go from summable digits to digits
aliroot -q -b "$ALICE_ROOT/ITS/AliITSSDigits2Digits.C"

# create reconstructed points for the ITS
aliroot -q -b "$ALICE_ROOT/ITS/AliITSDigits2RecPoints.C(0,$[$1-1])"

# prepare for tracking
aliroot -q -b "$ALICE_ROOT/ITS/AliITSTracksV1.C(0,$[$1-1])"

# do the tracking
aliroot -q -b "$ALICE_ROOT/ITS/AliITSTrackingV1.C(0,$[$1-1])"

# do the comparison
aliroot -q -b "$ALICE_ROOT/ITS/AliITSComparisonV1.C(0,$[$1-1])"

#
# after all of the above you can run AliITSPlotTracksV1.C macro under 
# aliroot to see plots of the efficiency and track parameter resolution
# or AliITSDisplayTracksV1.C to display found tracks on top of the
# ITS detailed geometry
