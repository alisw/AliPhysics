#!/bin/bash
CURDIR=`pwd`
cd $ALICE_ROOT
echo 'Making sure that TFluka is up to date...'
make all-TFluka
cd $CURDIR
# Remove the previous temporary directory
rm -rf tmp
# Make a new temporary directory and move to it
mkdir tmp
cd tmp
mkdir peg

# Link here some special Fluka files needed
ln -s $FLUPRO/xnloan.dat .
ln -s $FLUPRO/sigmapi.bin .
ln -s $FLUPRO/nuclear.bin .
ln -s $FLUPRO/neuxsc_72.bin neuxsc.bin
ln -s $FLUPRO/fluodt.dat .
ln -s $FLUPRO/elasct.bin .

# Copy the random seed
cp $FLUPRO/random.dat old.seed

# Give some meaningfull name to the output
ln -s fluka.out fort.11

# Link the pemf and input file for alice
ln -s $ALICE_ROOT/TFluka/input/alice.pemf .
#ln -s $ALICE_ROOT/TFluka/input/alice.inp .

#Link FlukaConfig.C as Config.C
ln -fs $ALICE_ROOT/TFluka/macro/FlukaConfig.C Config.C
echo 'Execute: gAlice->Init() OR gAlice->RunMC() at the ROOT prompt'
# Launch aliroot
aliroot
# Go back on exit
cd ..
