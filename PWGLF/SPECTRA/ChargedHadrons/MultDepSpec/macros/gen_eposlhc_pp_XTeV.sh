#!/bin/bash
set -e

if [[ $# -lt 2 ]]; then
    echo "You must specify at least where to write output to as first parameter and the center-of-mass energy in GeV as second parameter."
    exit -1
fi
beamErg=$(($2 / 2))

echo "running in `pwd`, writing HepMC to $1, pp collisions with center-of-mass energy $2 GeV"

# prepare environment
source /cvmfs/alice.cern.ch/etc/login.sh
eval $(alienv printenv CRMC::v1.5.4-3)

# force path to requested output
sed -e s#__CRMC_BASEDIR__#"$CRMC_ROOT"# ${ALICE_PHYSICS}/PWG/MCLEGO/CRMC/crmc_template.param > crmc.param

# run CRMC
crmc -o hepmc -p $beamErg -P-$beamErg -n 10000000 -m 0 -f $1
