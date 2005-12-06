#!/bin/bash
CURDIR=`pwd`
cd $ALICE_ROOT
echo 'Making sure that TFluka is up to date...'
make all-TFluka
cd $CURDIR

# Make working directory 
rm -rf fluka
mkdir -p fluka
cd fluka 

# Make a peg directory 
mkdir -p peg

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
# This is wrong:
#   ln -s $ALICE_ROOT/TFluka/input/FlukaVmc.pemf .
# Maybe 
cp $ALICE_ROOT/TFluka/input/alice.pemf FlukaVmc.pemf

#Link FlukaConfig.C as Config.C
cp $ALICE_ROOT/FMD/scripts/ConfigInner.C .
cp $ALICE_ROOT/.rootrc . 
# echo 'Execute: gAlice->Init() OR gAlice->RunMC() at the ROOT prompt'
# Launch aliroot
aliroot -l  # -b -q ../runIt.C > run.log 2>&1 

cd $CURDIR
# ____________________________________________________________________
#
# EOF
#
