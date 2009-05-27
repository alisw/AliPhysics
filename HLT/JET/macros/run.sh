#/bin/bash


if [ ! -d ./analysis ] ; then
    mkdir analysis
fi

ROOTFILES=`find . -maxdepth 1 -name "*.root" `


if [ -n "${ROOTFILES}" ] ; then
    rm ${ROOTFILES}
fi

aliroot -b -l -q 'HLTJetReconstruction.C(100,0,kTRUE)'





#valgrind --error-limit=no --leak-check=full --show-reachable=yes aliroot -b -l -q 'HLTJetReconstruction.C(10)' 2>&1 | tee valgrind.log
