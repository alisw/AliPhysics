#!/bin/bash

# -------------------------------------------
# Test VZERO reconstruction
# Author Jochen Thaeder <jochen@thaeder.de>
# -------------------------------------------

# N events
NEVENTS=20

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
aliroot -l -q -b ${ALICE_ROOT}/HLT/T0/macros/makeConfigurationObjectT0Reconstruction.C

# -- run chain for raw.root file
aliroot -l -q -b $ALICE_ROOT/HLT/T0/macros/HLTT0Test.C'("'${RAWPATH}'/raw.root","local://$ALICE_ROOT/OCDB",1,'${NEVENTS}')' 2>&1 | tee out.log

popd > /dev/null
