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
gGeoManager->Export("geometry.root");
.L $ALICE_ROOT/MUON/MUONCheckMisAligner.C+
MUONCheckMisAligner(0., 0.03, 0., 0.03, 0., 0.03, "FullMisAlignCDB");
.q
EOF

echo "Running simulation  ..."

aliroot -b  >& testSim.out << EOF 
// Uncoment following lines to run simulation with local residual mis-alignment
// (generated via MUONGenerateGeometryData.C macro)
AliCDBManager* man = AliCDBManager::Instance();
man->SetDefaultStorage("local://$ALICE_ROOT");
man->SetSpecificStorage("MUON/Align/Data","local://FullMisAlignCDB");
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
AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k5kG);
AliTracker::SetFieldMap(field, kFALSE);
AliReconstruction MuonRec("galice.root");
MuonRec.SetInput("$FULLPATH/");
MuonRec.SetRunTracking("MUON");
MuonRec.SetRunVertexFinder(kFALSE);
MuonRec.SetRunLocalReconstruction("MUON");
MuonRec.SetFillESD("MUON");
MuonRec.SetLoadAlignData("MUON")
MuonRec.Run();
.q
EOF

echo "Running alignment ..."

aliroot -b >& testAlign.out << EOF
.L $ALICE_ROOT/MUON/MUONAlignment.C+
MUONAlignment();
.q
EOF

echo "Finished"  
echo "... see results in test_align"

cd $CURDIR
