#!/bin/sh
# Remove the previous temporary directory
rm -rf tmp
# Make a new temporary directory and move to it
mkdir tmp
cd tmp

# Copy the random seed
cp $FLUPRO/data/random.dat old.seed

# Give some meaningfull name to the output
ln -s fluka.out fort.11
ln -s fluka.err fort.15
#Link FlukaConfig.C as Config.C
ln -fs $ALICE_ROOT/TFluka/macro/FlukaConfig.C Config.C
ln -fs $ALICE_ROOT/TFluka/input/coreFlukaVmc.inp coreFlukaVmc.inp
ln -fs $ALICE_ROOT/TFluka/macro/sim.C sim.C
echo 'Execute at the root prompt:'
echo '.x sim.C'

# Launch aliroot
aliroot
# Go back on exit
cd ..
