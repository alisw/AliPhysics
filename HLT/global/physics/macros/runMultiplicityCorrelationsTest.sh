#!/bin/bash

# -------------------------------------------
# Test MultiplictyCorrelations
# Author Jochen Thaeder <jochen@thaeder.de>
# -------------------------------------------

# N events
NEVENTS=60

# Path to ESD
#ESDPATH="/lustre/alice/jthaeder/data/compressionSGE/test2010"
ESDPATH="/lustre/alice/jthaeder/data/compressionSGE/test2011"
# -------------------------------------------

pushd $ESDPATH > /dev/null

#rm *.root  2> /dev/null
rm *.log   2> /dev/null
rm *.ps    2> /dev/null

if [ ! -d ./analysis ] ; then
    mkdir analysis
else
    rm ./analysis/*
fi

# -- Create config CDB object 
aliroot -l -q -b ${ALICE_ROOT}/HLT/global/physics/macros/makeConfigurationObjectMultiplicityCorrelations.C 2>&1 | tee configLog.log

# -- run chain for raw.root file
aliroot -l -q -b $ALICE_ROOT/HLT/global/physics/macros/HLTMultiplicityCorrelationsTest.C'("'${ESDPATH}'","local://$ALICE_ROOT/OCDB",'${NEVENTS}')' 2>&1 | tee out.log

popd > /dev/null
