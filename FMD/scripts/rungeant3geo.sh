#!/bin/bash
CURDIR=`pwd`
cd $ALICE_ROOT
echo 'Making sure that AliROOT is up to date...'
make
cd $CURDIR

# Make working directory 
rm -rf geant321
mkdir -p geant321
cd geant321


#Link FlukaConfig.C as Config.C
cp $ALICE_ROOT/FMD/Config.C .
cp $ALICE_ROOT/.rootrc . 
# echo 'Execute: gAlice->Init() OR gAlice->RunMC() at the ROOT prompt'
# Launch aliroot
aliroot -l -b -q ../runIt.C > run.log 2>&1 

cd $CURDIR
