#!/bin/sh
# $Id$
# A. De Falco, H. Woehri, INFN Cagliari, April 2007

NEVENTS=1
CURDIR=`pwd`

declare -i time
time=`date +%s` 
SEED=$time 
SEED2=$RANDOM 
OUTDIR=fastOut_$SEED$SEED2

rm -fr $OUTDIR
mkdir $OUTDIR
cp $ALICE_ROOT/MUON/.rootrc $ALICE_ROOT/MUON/rootlogon.C $OUTDIR
cd $OUTDIR
#echo current dircetory is $PWD
  
###############################
#  Event Generation           #
###############################
echo 'Performing the fast generation now ...'
aliroot -b > gen.out 2> gen.err << EOF  
gSystem->Load("libFASTSIM");
.L $ALICE_ROOT/MUON/fastMUONGen.C+
fastMUONGen($NEVENTS, "galice.root", 2);
.q
EOF

#####################################
#  Event Simulation (Fast Tracking) #
#####################################
echo 'Performing the fast reconstruction now ...'
aliroot -b > sim.out 2> sim.err << EOF  
.L $ALICE_ROOT/MUON/fastMUONSim.C+
fastMUONSim();
.q
EOF
