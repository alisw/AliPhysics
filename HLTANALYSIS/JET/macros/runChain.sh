#/bin/bash

NEVENTS=200
Pt_BIN=kPythia6Jets104_125

if [ ! -d ./analysis ] ; then
    mkdir analysis
fi

ROOTFILES=`find . -maxdepth 1 -name "*.root" `

if [ -n "${ROOTFILES}" ] ; then
    rm ${ROOTFILES}
fi
 

aliroot -b -l -q HLTJetReconstruction.C'('${NEVENTS}', kFALSE, '${Pt_BIN}')' 2>&1 | tee ChainLog.log

ROOTFILES=`find . -maxdepth 1 -name "*.root" `

if [ -n "${ROOTFILES}" ] ; then
    rm ${ROOTFILES}
fi

aliroot -l -b -q readJets.C'("./analysis/EOR_analyze_'${NEVENTS}'_'${Pt_BIN}'.root")'

# valgrind --error-limit=no --leak-check=full --show-reachable=yes aliroot -b -l -q 'HLTJetReconstruction.C(10)' 
