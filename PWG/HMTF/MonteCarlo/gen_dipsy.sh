#!/bin/bash
set -e

if [[ $# -lt 1 ]]; then
    echo "You must specify at least where to write output to as first parameter."
    exit -1
fi

echo "running in `pwd`, writing HepMC to $1"

# force path to requested output
sed -e "s/^\([[:space:]]*set HepMCFile:Filename\) .*$/\1 $1/" ${ALICE_PHYSICS}/PWG/HMTF/MonteCarlo/DIPSYpp_HepMC.in > thepeg.in
# setup DIPSY
setupThePEG -r ${THEPEG_BASEDIR}/lib/ThePEG/ThePEGDefaults.rpo -I ${THEPEG_BASEDIR}/share/Ariadne thepeg.in > setup.log 2>&1
# run DIPSY
runThePEG DIPSYpp.run -N 2147483647
