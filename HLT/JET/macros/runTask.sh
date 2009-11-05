#/bin/bash

if [ ! -d ./analysis ] ; then
    mkdir analysis
fi

ROOTFILES=`find . -maxdepth 1 -name "*.root" `

if [ -n "${ROOTFILES}" ] ; then
    rm ${ROOTFILES}
fi

aliroot -b -l -q './tasks/JetAnalysisManagerHLT.C' 2>&1 | tee TaskLog.log

ROOTFILES=`find . -maxdepth 1 -name "*.root" `

#if [ -n "${ROOTFILES}" ] ; then
#    rm ${ROOTFILES}
#fi

#aliroot -l 'readJets.C("./analysis/EOR_analyze_100_kPythia6Jets104_125.root")'


#valgrind --error-limit=no --leak-check=full --show-reachable=yes aliroot -b -l -q 'HLTJetReconstruction.C(10)' 
