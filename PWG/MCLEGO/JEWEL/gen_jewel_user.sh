#!/bin/bash
set -e

if [[ $# -lt 3 ]]; then
    echo "You must specify exactly 3 parameters."
    exit -1
fi

echo "running in `pwd`, writing HepMC to $1"

# prepare environment
source /cvmfs/alice.cern.ch/etc/login.sh
eval $(alienv printenv JEWEL::v2.0.2-2)

# take configuration from arguments
echo -e "$2" >> params-simple.dat
echo -e "$3" >> medium-params.dat

# Add NJOB to time-based random seed
nrandom=`date '+%d%H%M%S'`
sed -i "1s/^/NJOB $nrandom\n/" ./params-simple.dat

# run JEWEL
jewel-2.0.2-simple params-simple.dat
