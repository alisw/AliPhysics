#!/bin/sh
# $Id$
#
# By J. Castillo


CURDIR=`pwd`
OUTDIR=test_align

rm -fr $OUTDIR
mkdir $OUTDIR
cp $ALICE_ROOT/MUON/.rootrc $ALICE_ROOT/MUON/rootlogon.C $OUTDIR
cd $OUTDIR

FULLPATH="$CURDIR/$OUTDIR"
NEVENTS=1000
SEED=1234567

echo "Generating misalignment ..."

aliroot -b >& testMisalign.out << EOF
gAlice->Init("$ALICE_ROOT/MUON/Config.C");
.x $ALICE_ROOT/MUON/MUONCheckMisAligner.C(0., 0.03, 0., 0.03, 0., 0.03, "FullMisAlignCDB");
.q
EOF

cp -r $ALICE_ROOT/MUON/Calib FullMisAlignCDB/MUON/

echo "Running simulation  ..."

aliroot -b  >& testSim.out << EOF 
// Uncoment following lines to run simulation with local residual mis-alignment
// (generated via MUONGenerateGeometryData.C macro)
AliCDBManager* man = AliCDBManager::Instance();
man->SetDefaultStorage("local://$ALICE_ROOT");
man->SetSpecificStorage("MUON","local://FullMisAlignCDB");
gRandom->SetSeed($SEED);
AliSimulation MuonSim("$ALICE_ROOT/MUON/Config.C");
MuonSim.SetMakeTrigger("MUON");
MuonSim.SetWriteRawData("MUON");
MuonSim.Run($NEVENTS);
.q
EOF

echo "Running reconstruction  ..."

aliroot -b >& testReco.out << EOF
gRandom->SetSeed($SEED);
AliReconstruction MuonRec("galice.root");
MuonRec.SetInput("$FULLPATH/");
MuonRec.SetRunTracking("");
MuonRec.SetRunVertexFinder(kFALSE);
MuonRec.SetRunLocalReconstruction("MUON");
MuonRec.SetFillESD("MUON");
MuonRec.SetLoadAlignData("MUON")
MuonRec.Run();
.q
EOF

echo "Running alignment ..."

aliroot -b >& testAlign.out << EOF
.x $ALICE_ROOT/MUON/MUONAlignment.C
.q
EOF

echo "Finished"  
echo "... see results in test_align"

cd $CURDIR
