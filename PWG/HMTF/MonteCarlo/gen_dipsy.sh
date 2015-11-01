#!/bin/bash
set -e

if [[ $# -lt 1 ]]; then
    echo "You must specify at least where to write output to as first parameter."
    exit -1
fi

echo "running in `pwd`, writing HepMC to $1"

# prepare environment
source /cvmfs/alice.cern.ch/etc/login.sh
eval $(alienv printenv ThePEG::v2015-08-11-3)

# force path to requested output
sed -e "s/^\([[:space:]]*set HepMCFile:Filename\) .*$/\1 $1/" ${ALICE_PHYSICS}/PWG/HMTF/MonteCarlo/DIPSYpp_HepMC.in > thepeg.in

# setup DIPSY
setupThePEG -r ${THEPEG_ROOT}/lib/ThePEG/ThePEGDefaults.rpo -I ${THEPEG_ROOT}/share/Ariadne thepeg.in > setup.log 2>&1

# run DIPSY
PYTHIA8DATA=$PYTHIA_ROOT/share/Pythia8/xmldoc runThePEG DIPSYpp.run -N 2147483647 --seed ${ALIEN_PROC_ID:=0}
