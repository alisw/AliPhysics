#!/bin/sh
# Remove the previous temporary directory
rm -rf tmp
# Make a new temporary directory and move to it
mkdir tmp
cd tmp

# Link here some special Fluka files needed
ln -s $FLUPRO/data/neuxsc.bin  .
ln -s $FLUPRO/data/elasct.bin  .
ln -s $FLUPRO/data/gxsect.bin  .
ln -s $FLUPRO/data/neuxsc.bin  .
ln -s $FLUPRO/data/nuclear.bin .
ln -s $FLUPRO/data/sigmapi.bin .
ln -s $FLUPRO/data/brems_fin.bin .
ln -s $FLUPRO/data/cohff.bin .
ln -s $FLUPRO/data/fluodt.dat  .
ln -s $FLUPRO/data/random.dat  .
# Copy the random seed
cp $FLUPRO/random.dat old.seed

# Give some meaningfull name to the output
ln -s fluka.out fort.11
ln -s fluka.err fort.15
#Link FlukaConfig.C as Config.C
ln -fs $ALICE_ROOT/TFluka/macro/FlukaConfig.C Config.C
ln -fs $ALICE_ROOT/TFluka/input/coreFlukaVmc.inp coreFlukaVmc.inp
ln -fs $ALICE_ROOT/TFluka/macro/sim.C sim.C
echo 'Execute at the root prompt:'
echo 'AliSimulation sim'
echo 'sim.Run()'

# Launch aliroot
aliroot
# Go back on exit
cd ..
