#!/bin/bash
#./runExtractGenMCtruth.sh "finalCuts/MC_pp/7TeV/LHC10f6a/corrected/finalisedSplines/analytical/Jets/nclCut/bhess_PID_Jets.root" "" 0 -2 -2

aliroot 'extractGenMCtruth.C+("'$1'", "'$2'",  '$3', '$4', '$5', 5, 10)' -l -b -q

aliroot 'extractGenMCtruth.C+("'$1'", "'$2'",  '$3', '$4', '$5', 10, 15)' -l -b -q

aliroot 'extractGenMCtruth.C+("'$1'", "'$2'",  '$3', '$4', '$5', 15, 20)' -l -b -q

aliroot 'extractGenMCtruth.C+("'$1'", "'$2'",  '$3', '$4', '$5', 20, 40)' -l -b -q

aliroot 'extractGenMCtruth.C+("'$1'", "'$2'",  '$3', '$4', '$5', 40, 80)' -l -b -q

aliroot 'extractGenMCtruth.C+("'$1'", "'$2'",  '$3', '$4', '$5', 10, 80)' -l -b -q

aliroot 'extractGenMCtruth.C+("'$1'", "'$2'",  '$3', '$4', '$5', 20, 30)' -l -b -q

aliroot 'extractGenMCtruth.C+("'$1'", "'$2'",  '$3', '$4', '$5', 30, 40)' -l -b -q

aliroot 'extractGenMCtruth.C+("'$1'", "'$2'",  '$3', '$4', '$5', 10, 40)' -l -b -q
