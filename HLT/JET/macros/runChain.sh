#/bin/bash

if [ ! -d ./analysis ] ; then
    mkdir analysis
fi

ROOTFILES=`find . -maxdepth 1 -name "*.root" `

if [ -n "${ROOTFILES}" ] ; then
    rm ${ROOTFILES}
fi

aliroot -b -l -q 'HLTJetReconstruction.C(20,kTRUE)' 2>&1 | tee ChainLog.log

ROOTFILES=`find . -maxdepth 1 -name "*.root" `

if [ -n "${ROOTFILES}" ] ; then
    rm ${ROOTFILES}
fi

aliroot -l 'readJets.C("./analysis/EOR_analyze_20_kPythia6Jets104_125.root")'

#valgrind --error-limit=no --leak-check=full --show-reachable=yes aliroot -b -l -q 'HLTJetReconstruction.C(10)' 
