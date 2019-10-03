#!/bin/bash

### configuration
FIFOPATH=/tmp/DIPSYpp.hepmc
NEV=10000

### execution
mkfifo ${FIFOPATH}
setupThePEG -r ${THEPEG_BASEDIR}/lib/ThePEG/ThePEGDefaults.rpo -I ${THEPEG_BASEDIR}/share/Ariadne DIPSYpp_HepMC.in
runThePEG DIPSYpp.run -N${NEV} --tics &
./run.sh -t kHepMC -i ${FIFOPATH} -n ${NEV} -r both
rm ${FIFOPATH}
