#!/bin/sh
DeleteOldFiles
# do simulation, digitization, reconstruction and tracking in TPC
aliroot -q AliTPCtracks.C
# do digitization in ITS
aliroot -q ITSHitsToDigits.C
# do reconstruction in ITS
aliroot -q ITSDigitsToClusters.C
# do tracking in ITS
aliroot -q ITStracks.C
aliroot -q ITStracking.C
