#!/bin/bash
set -e

if [[ $# -lt 1 ]]; then
    echo "You must specify at least where to write output to as first parameter."
    exit -1
fi

echo "running in `pwd`, writing HepMC to $1"

# prepare environment
source /cvmfs/alice.cern.ch/etc/login.sh
eval $(alienv printenv AMPT::v1.26t7-v2.26t7-1)

# run generator

# generate random seed
nrandom=`date '+%d%H%M%S'`
echo $nrandom > nseed_runtime

if [ -d ana/ ];
then
    rm -rf ana/
fi
mkdir ana
cp ${ALICE_PHYSICS}/PWG/MCLEGO/AMPT/input.ampt .
cp input.ampt ana/

mkfifo ana/ampt.dat

ampt < nseed_runtime &
parser ana/ampt.dat $1

rm -rf ana/
