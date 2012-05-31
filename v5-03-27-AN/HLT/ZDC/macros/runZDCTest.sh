#!/bin/bash

# -------------------------------------------
# Test ZDC reconstruction
# Author Chiara Oppedisano <Chiara.Oppedisano@to.infn.it>
#        Jochen Thaeder <jochen@thaeder.de>
# -------------------------------------------

# N events
NEVENTS=1

# Path to raw.root
RAWPATH="/opt/HLT/aliroot/AliRoot_HEAD_2010-09-01/test/ppbench/recraw"

# -------------------------------------------

pushd $RAWPATH > /dev/null

rm *.root  2> /dev/null
rm *.log   2> /dev/null
rm *.ps    2> /dev/null

ln -s ../raw.root

if [ ! -d ./analysis ] ; then
    mkdir analysis
else
    rm ./analysis/*
fi

# -- Create config CDB object 
#aliroot -l -q -b ${ALICE_ROOT}/HLT/ZDC/macros/makeConfigurationObjectZDCReconstruction.C

# -- run chain for raw.root file
valgrind --tool=callgrind aliroot -l -q -b $ALICE_ROOT/HLT/ZDC/macros/HLTZDCTest.C'("'${RAWPATH}'/raw.root","local://$ALICE_ROOT/OCDB",1,'${NEVENTS}')' 2>&1 | tee out.log

popd > /dev/null