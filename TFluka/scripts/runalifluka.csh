#####################################################################
# runalifluka.csh: Script File that sets up the needed input/output
# files for TFluka
# Author: I. Gonzalez
# Date: 11/2002
#####################################################################
# Remove the previous temporary directory
rm -rf tmp
# Make a new temporary directory and move to it
mkdir tmp
cd tmp

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

# Launch aliroot
aliroot -b <<EOF > run.out
gAlice->Init();
gAlice->Run();
EOF

# Go back on exit
cd ..
