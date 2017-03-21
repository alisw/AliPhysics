#!/bin/bash

# -------------------------------------------
# Test T0 reconstruction
# Author Jochen Thaeder <jochen@thaeder.de>
# -------------------------------------------

# N events
NEVENTS=1000

# Path to raw.root
#RAWPATH="/opt/HLT/aliroot/AliRoot_HEAD_2010-09-01/test/ppbench/recraw"
RAWPATH="/tzero/alla/alice/AliRootHLT/HLT/T0/macros"

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
#aliroot -l -q -b $ALICE_ROOT/HLT/T0/macros/HLTT0Test.C'("'${RAWPATH}'/raw.root","raw://",1,'${NEVENTS}')' 2>&1 | tee out.log
aliroot -l -q -b $ALICE_ROOT/HLT/T0/macros/T0CalibrationProcessor.C'("alien:///alice/data/2012/LHC12h/000192729/raw/12000192729062.55.root","raw://",1,'${NEVENTS}')' 2>&1 | tee out.log
#aliroot -l -q -b $ALICE_ROOT/HLT/T0/macros/HLTT0Test.C'("alien:///alice/data/2012/LHC12h/000192729/raw/12000192729062.55.root","raw://",1,'${NEVENTS}')' 2>&1 | tee out.log

popd > /dev/null
