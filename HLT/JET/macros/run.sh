#/bin/bash


N_EVENTS=100

if [ ! -d ./analysis ] ; then
    mkdir analysis
fi

ROOTFILES=`find . -maxdepth 1 -name "*.root" `

if [ -n "${ROOTFILES}" ] ; then
    rm ${ROOTFILES}
fi


aliroot -b -l -q 'HLTJetReconstruction.C(10000,0,kTRUE)' 2>&1 | tee log.log

ROOTFILES=`find . -maxdepth 1 -name "*.root" `

if [ -n "${ROOTFILES}" ] ; then
    rm ${ROOTFILES}
fi

aliroot -l 'readJets.C("./analysis/EOR_analyze_10000_kPythia6Jets104_125.root")'


#valgrind --error-limit=no --leak-check=full --show-reachable=yes aliroot -b -l -q 'HLTJetReconstruction.C(10)' 
