#!/bin/bash
set -e

if [[ $# -lt 1 ]]; then
    echo "You must specify at least where to write output to as first parameter."
    exit -1
fi

echo "running in `pwd`, writing HepMC to $1"

# prepare environment
# source /cvmfs/alice.cern.ch/etc/login.sh

if [ -d events/ ];
then
    rm -rf events/
fi
mkdir events
mkdir events/lhyquid2dbi-RHICAuAu200c2030Ti455ti025Tf145

cp ${ALICE_PHYSICS}/PWG/MCLEGO/Therminator2/events.ini .
cp -r $THERMINATOR2_ROOT/fomodel .
cp -r $THERMINATOR2_ROOT/share .
cp -r $THERMINATOR2_ROOT/macro .

mkfifo events/lhyquid2dbi-RHICAuAu200c2030Ti455ti025Tf145/event.txt

therm2_events events.ini 0 &
therm2_parser events/lhyquid2dbi-RHICAuAu200c2030Ti455ti025Tf145/event.txt $1 

