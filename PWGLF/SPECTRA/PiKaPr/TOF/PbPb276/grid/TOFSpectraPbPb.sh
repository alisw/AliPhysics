#! /bin/bash

tar zxvf TOFSpectraPbPb.tgz

if [ "$1" == "data" ]; then
    echo "*** Running analysis - Data ***"
    aliroot -b -q "SteerAnalysisTaskTOFSpectraPbPb.C(\"wn.xml\", kFALSE, kFALSE, kTRUE)" 2>&1 | tee TOFSpectraPbPb.log
elif [ "$1" == "mc" ]; then
    echo "*** Running analysis - MC ***"
    aliroot -b -q "SteerAnalysisTaskTOFSpectraPbPb.C(\"wn.xml\", kTRUE, kFALSE, kTRUE)" 2>&1 | tee TOFSpectraPbPb.log
elif [ "$1" == "mctune" ]; then
    echo "*** Running analysis - MCTune ***"
    aliroot -b -q "SteerAnalysisTaskTOFSpectraPbPb.C(\"wn.xml\", kTRUE, kTRUE, kTRUE)" 2>&1 | tee TOFSpectraPbPb.log
else
    echo "*** Unkwnown analysis ***"
    aliroot -b -q "SteerAnalysisTaskTOFSpectraPbPb.C(\"wn.xml\", kFALSE, kFALSE, kTRUE)" 2>&1 | tee TOFSpectraPbPb.log
fi
