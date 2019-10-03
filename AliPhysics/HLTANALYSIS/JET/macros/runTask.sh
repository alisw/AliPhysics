#/bin/bash

if [ "$1" == "MC" ] ; then 
    MODE=MC 
elif [ "$1" == "ESD" ] ; then 
    MODE=ESD
else    
    MODE=MC
fi

ROOTFILES=`find . -maxdepth 1 -name "*.root" `

if [ -n "${ROOTFILES}" ] ; then
    rm ${ROOTFILES}
fi

if [ "$MODE" == "ESD" ] ; then
    echo " -= Load : JetAnalysisManagerHLT.C =- "
    aliroot -b -l -q './tasks/JetAnalysisManagerHLT.C' 2>&1 | tee TaskLog.log
else
    echo " -= Load : JetAnalysisManagerHLTMC.C =- "
    aliroot -b -l -q './tasks/JetAnalysisManagerHLTMC.C' 2>&1 | tee TaskLog.log
fi

ROOTFILES=`find . -maxdepth 1 -name "*.root" `

#if [ -n "${ROOTFILES}" ] ; then
#    rm ${ROOTFILES}
#fi

#aliroot -l 'readJets.C("./analysis/EOR_analyze_100_kPythia6Jets104_125.root")'


#valgrind --error-limit=no --leak-check=full --show-reachable=yes aliroot -b -l -q 'HLTJetReconstruction.C(10)' 
