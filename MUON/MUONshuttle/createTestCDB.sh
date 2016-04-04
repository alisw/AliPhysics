#!/bin/sh

# Create a set of links under $ALICE_ROOT/../src/SHUTTLE/TestShuttle/TestCDB
# to be able to use the TestMUONPreprocessor.C macro

rm -rf $ALICE_ROOT/../src/SHUTTLE/TestShuttle/TestCDB/*

cd $ALICE_ROOT/../src/SHUTTLE/TestShuttle/TestCDB
mkdir -p MUON/Calib/MappingData
cd MUON/Calib/MappingData/
ln -si $ALICE_ROOT/OCDB/MUON/Calib/MappingData/* .
cd $ALICE_ROOT/../src/SHUTTLE/TestShuttle/TestCDB
mkdir -p MUON/Calib/Config
cd MUON/Calib/Config
ln -si $ALICE_ROOT/OCDB/MUON/Calib/Config/* .

cd $ALICE_ROOT/../src/SHUTTLE/TestShuttle/TestCDB
mkdir -p MUON/Align
cd MUON/Align
ln -si $ALICE_ROOT/OCDB/MUON/Align/Baseline .
