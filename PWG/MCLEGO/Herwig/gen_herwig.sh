#!/bin/bash
set -e

if [[ $# -lt 1 ]]; then
    echo "You must specify at least where to write output to as first parameter."
    exit -1
fi

echo "running in `pwd`, writing HepMC to $1"

# prepare environment
source /cvmfs/alice.cern.ch/etc/login.sh
eval $(alienv printenv Herwig::v7.0.4-alice1-1)

# run Herwig
sed -e s#output.hepmc#"$1"# ${ALICE_PHYSICS}/PWG/MCLEGO/Herwig/LHC-MB.in > Matchbox.in
Herwig read --repo=$HERWIG_ROOT/share/Herwig/HerwigDefaults.rpo Matchbox.in
