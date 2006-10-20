#!/bin/sh
# $Id$

CURDIR=`pwd`
OUTDIR=test_out

rm -fr $OUTDIR
mkdir $OUTDIR
cp $ALICE_ROOT/MUON/.rootrc $ALICE_ROOT/MUON/rootlogon.C $OUTDIR
cd $OUTDIR

FULLPATH="$CURDIR/$OUTDIR"
NEVENTS=100
SEED=1234567

echo "Running simulation  ..."

aliroot -b  >& testSim.out << EOF 
// Uncoment following lines to run simulation with local residual mis-alignment
// (generated via MUONGenerateGeometryData.C macro)
// AliCDBManager* man = AliCDBManager::Instance();
// man->SetDefaultStorage("local://$ALICE_ROOT");
// man->SetSpecificStorage("MUON","local://$ALICE_ROOT/MUON/ResMisAlignCDB");
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
AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k4kG);
AliTracker::SetFieldMap(field, kFALSE);
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

echo "Running Trigger efficiency  ..."
aliroot -b >& testTriggerResults.out << EOF
.L $ALICE_ROOT/MUON/MUONTriggerEfficiency.C+
MUONTriggerEfficiency();
.q
EOF

echo "Running efficiency  ..."

aliroot -b >& testResults.out << EOF
.L $ALICE_ROOT/MUON/MUONefficiency.C+
// no argument assumes Upsilon but MUONefficiency(443) works on Jpsi
MUONefficiency();
.q
EOF

if [ "$NEVENTS" -le 20 ]; then

echo "Running dumps ..."

aliroot -b << EOF
.L $ALICE_ROOT/MUON/MUONCheck.C+
MUONdigits(); > check.digits
MUONrecpoints(); > check.recpoints
MUONrectracks(); > check.rectracks
MUONrectrigger(); > check.rectrigger
.q
EOF

fi

echo "Finished"  
echo "... see results in test_out"

cd $CURDIR
