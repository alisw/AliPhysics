#!/bin/sh
# delete eventual old files from the last run
DeleteOldFiles
# run the hit generation
aliroot -q -b grun.C
# digitize TPC
aliroot -q -b AliTPCHits2Digits.C
# do tracking in TPC
aliroot -q -b AliTPCtracknew.C
# prepare TPC tracks for matching with the ITS
aliroot -q -b TPCtracks.C
# digitize ITS
aliroot -q -b ITSHitsToDigits.C
# create reconstruct point for the ITS
aliroot -q -b ITSDigitsToClusters.C
# prepare for tracking
aliroot -q -b ITStracks.C
# do the tracking
aliroot -q -b ITStracking.C
#
# after all of the above you can run ITSPlotTracks.C macro under aliroot to
# see plots of the efficiency and track parameter resolution
