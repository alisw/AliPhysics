#!/bin/sh
# $Id$
#
# By J. Castillo

CURDIR=`pwd`
OUTDIR=test_align
SIMDIR="generated"

rm -fr $OUTDIR
mkdir $OUTDIR
cp $ALICE_ROOT/MUON/.rootrc $ALICE_ROOT/MUON/rootlogon.C $OUTDIR
cd $OUTDIR

FULLPATH="$CURDIR/$OUTDIR"
NEVENTS=1000
SEED=1234567

echo "Generating misalignment ..."

aliroot -b >& testMisalign.out << EOF
AliMpCDB::LoadMpSegmentation2();
TString configFileName = "$ALICE_ROOT/MUON/Config.C";
gROOT->LoadMacro(configFileName.Data());
gInterpreter->ProcessLine(gAlice->GetConfigFunction());
gAlice->GetMCApp()->Init();
gGeoManager->Export("geometry.root");
.L $ALICE_ROOT/MUON/MUONCheckMisAligner.C+
MUONCheckMisAligner(0., 0.03, 0., 0.03, 0., 0.03, "FullMisAlignCDB");
.q
EOF

echo "Running simulation  ..."

aliroot -b  >& testSim.out << EOF 
// Uncoment following lines to run simulation with local full mis-alignment
// (generated via MUONGenerateGeometryData.C macro)
AliCDBManager* man = AliCDBManager::Instance();
man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
man->SetSpecificStorage("MUON/Align/Data","local://FullMisAlignCDB");
AliSimulation MuonSim("$ALICE_ROOT/MUON/Config.C");
MuonSim.SetSeed($SEED);
MuonSim.SetMakeTrigger("MUON");
MuonSim.SetWriteRawData("MUON","raw.root",kTRUE);
MuonSim.SetMakeDigits("MUON");
MuonSim.SetMakeSDigits("MUON");
MuonSim.SetMakeDigitsFromHits("");
MuonSim.Run($NEVENTS);
.q
EOF

echo "Moving generated files to $SIMDIR"
mkdir $SIMDIR
mv *QA*.root *.log $SIMDIR
mv MUON*.root Kinematics*.root galice.root TrackRefs*.root $SIMDIR

echo "Running reconstruction  ..."

aliroot -b >& testReco.out << EOF
// Uncoment following lines to run reconstruction with local full mis-alignment
// (generated via MUONGenerateGeometryData.C macro)
//AliCDBManager* man = AliCDBManager::Instance();
//man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
//man->SetSpecificStorage("MUON/Align/Data","local://$ALICE_ROOT/OCDB/MUON/FullMisAlignCDB");
gRandom->SetSeed($SEED);
AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k5kG);
AliTracker::SetFieldMap(field, kFALSE);
AliReconstruction MuonRec("galice.root");
MuonRec.SetInput("$FULLPATH/raw.root");
MuonRec.SetRunVertexFinder(kFALSE);
MuonRec.SetRunLocalReconstruction("MUON");
MuonRec.SetRunTracking("MUON");
MuonRec.SetFillESD("");
MuonRec.SetLoadAlignData("MUON")
MuonRec.SetNumberOfEventsPerFile(1000);
AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::GetLowFluxParam();
//  muonRecoParam->SaveFullClusterInESD(kTRUE,100.);
muonRecoParam->CombineClusterTrackReco(kFALSE);
//muonRecoParam->SetClusteringMode("PEAKFIT");
//muonRecoParam->SetClusteringMode("PEAKCOG");
muonRecoParam->Print("FULL");

AliRecoParam::Instance()->RegisterRecoParam(muonRecoParam);
MuonRec.Run();
.q
EOF

echo "Running alignment ..."

aliroot -b >& testAlign.out << EOF
.L $ALICE_ROOT/MUON/MUONAlignment.C+
AliMpCDB::LoadMpSegmentation2();
MUONAlignment();
.q
EOF

echo "Finished"  
echo "... see results in " $OUTDIR

cd $CURDIR
