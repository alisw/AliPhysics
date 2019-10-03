#!/bin/bash
set -e

if [[ $# -lt 1 ]]; then
    echo "You must specify at least where to write output to as first parameter."
    exit -1
fi

echo "running in `pwd`, writing HepMC to $1"

# prepare environment
source /cvmfs/alice.cern.ch/etc/login.sh
eval $(alienv printenv JEWEL::v2.0.2-2)

# prepare configuration
sed -e "s/^\([[:space:]]*HEPMCFILE\) .*$/\1 $1/" ${ALICE_PHYSICS}/PWG/MCLEGO/JEWEL/params-simple.dat > params-simple.dat
cp ${ALICE_PHYSICS}/PWG/MCLEGO/JEWEL/params.medium-simple.dat .

# run JEWEL
jewel-2.0.2-simple params-simple.dat
