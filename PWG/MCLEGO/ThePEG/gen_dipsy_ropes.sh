#!/bin/bash
set -e

if [[ $# -lt 1 ]]; then
    echo "You must specify at least where to write output to as first parameter."
    exit -1
fi

echo "running in `pwd`, writing HepMC to $1"

# prepare environment
source /cvmfs/alice.cern.ch/etc/login.sh
eval $(alienv printenv ThePEG::v2015-08-11-4)

# force path to requested output
sed -e "s/^\([[:space:]]*set HepMCFile:Filename\) .*$/\1 $1/" ${ALICE_PHYSICS}/PWG/MCLEGO/ThePEG/Rope.in > thepeg.in
cp ${ALICE_PHYSICS}/PWG/MCLEGO/ThePEG/LHCpp.in .
cp ${ALICE_PHYSICS}/PWG/MCLEGO/ThePEG/Tune31.in .

# setup DIPSY
setupThePEG -r ${THEPEG_ROOT}/lib/ThePEG/ThePEGDefaults.rpo -I ${THEPEG_ROOT}/share/Ariadne thepeg.in > setup.log 2>&1

# run DIPSY
runThePEG Rope.run -N 2147483647 --seed ${ALIEN_PROC_ID:=0}
