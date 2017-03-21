#!/bin/bash

# -------------------------------------------
# Test TZERO reconstruction
# Author Jochen Thaeder <jochen@thaeder.de>
# -------------------------------------------

# N events
NEVENTS=10
echo $NEVENTS
# Path to raw.root
#RAWPATH="/opt/HLT/aliroot/AliRoot_HEAD_2010-09-01/test/ppbench/recraw"
RAWPATH="/tzero/alla/alice/AliRootHLT/HLT/T0/macros"
echo $RAWPATH
# -------------------------------------------

#pushd $RAWPATH > /dev/null

#rm *.root  2> /dev/null
rm galice.root 
#rm *.log   2> /dev/null
#rm *.ps    2> /dev/null

#ln -s ../raw.root

#if [ ! -d ./analysis ] ; then
#    mkdir analysis
#else
#    rm ./analysis/*
#fi
# -- Create config CDB object 
aliroot -q  ${ALICE_ROOT}/HLT/T0/macros/makeConfigurationObjectTZEROReconstruction.C
# -- run chain for raw.root file
#aliroot -l -q -b $ALICE_ROOT/HLT/T0/macros/HLTTZEROTest.C'("'${RAWPATH}'/raw.root","raw://",1,'${NEVENTS}')' 2>&1 | tee out.log
aliroot -l -q -b $ALICE_ROOT/HLT/T0/macros/TZEROCalibrationProcessor.C'("/home/alla/alice/AliRootHLT/HLT/T0/macros/10000128483044.40.root","raw://",1,'${NEVENTS}')' 2>&1 | tee out.log
#aliroot -l -q -b $ALICE_ROOT/HLT/T0/macros/HLTTZEROCalibrationProcessor.C'("alien:///alice/data/2012/LHC12h/000192729/raw/12000192729062.55.root","raw://",1,'${NEVENTS}')' 2>&1 | tee out.log

#popd > /dev/null
