#!/bin/bash

set -e

if [[ $# -lt 1 ]]; then
    echo "You must specify at least where to write output to as first parameter."
    exit -1
fi

echo "running in `pwd`, writing HepMC to $1"

# prepare environment
source /cvmfs/alice.cern.ch/etc/login.sh
eval $(alienv printenv Therminator2::v2.0.3-1)

if [ -d events/ ];
then
    rm -rf events/
fi
mkdir events

cp ${ALICE_PHYSICS}/PWG/MCLEGO/Therminator2/events.ini .
cp ${ALICE_PHYSICS}/PWG/MCLEGO/Therminator2/apply_custom_config.py .
cp -r $THERMINATOR2_ROOT/fomodel .
cp -r $THERMINATOR2_ROOT/share .
cp -r $THERMINATOR2_ROOT/macro .

mkfifo events/event.txt

# set the key=value pairs in the configuration files
# also downloads and sets the proper xml file (if hydro model is selected and custom xml path set)
python apply_custom_config.py

# start the generator and the parser
therm2_parser events/event.txt $1 &
therm2_events events.ini 0
 

