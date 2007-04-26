#!/bin/sh
# $Id$

CURDIR=`pwd`
OUTDIR=test_out

rm -fr $OUTDIR
mkdir $OUTDIR
cp $ALICE_ROOT/MUON/.rootrc $ALICE_ROOT/MUON/rootlogon.C $OUTDIR
cd $OUTDIR

RUN=0
FULLPATH="$CURDIR/$OUTDIR"
NEVENTS=100
SEED=1234567

echo "Running simulation  ..."

aliroot -b  >& testSim.out << EOF 
// Uncoment following lines to run simulation with local residual mis-alignment
// (generated via MUONGenerateGeometryData.C macro)
// AliCDBManager* man = AliCDBManager::Instance();
// man->SetDefaultStorage("local://$ALICE_ROOT");
// man->SetSpecificStorage("MUON/Align/Data","local://$ALICE_ROOT/MUON/ResMisAlignCDB");
gRandom->SetSeed($SEED);
AliCDBManager::Instance()->SetRun($RUN);
AliSimulation MuonSim("$ALICE_ROOT/MUON/Config.C");
MuonSim.SetMakeTrigger("MUON");
MuonSim.SetWriteRawData("MUON");
MuonSim.Run($NEVENTS);
.q
EOF

echo "Removing Digits files ..."
mkdir MUON.Digits
mv MUON.Digits*.root MUON.Digits/ 

echo "Running reconstruction  ..."

aliroot -b >& testReco.out << EOF
AliCDBManager::Instance()->SetRun($RUN);
gRandom->SetSeed($SEED);
AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k5kG);
AliTracker::SetFieldMap(field, kFALSE);
AliReconstruction MuonRec("galice.root");
MuonRec.SetInput("$FULLPATH/");
MuonRec.SetRunVertexFinder(kFALSE);
MuonRec.SetRunLocalReconstruction("MUON");
MuonRec.SetRunTracking("MUON");
MuonRec.SetFillESD("MUON");
MuonRec.SetLoadAlignData("MUON");
// Uncoment following line to run reconstruction with the orginal tracking method
// instead of the kalman one (default)
//MuonRec.SetOption("MUON","Original");
// Use the following to change clustering method
//MuonRec.SetOption("MUON","MLEM"); // new scheme AZs clustering
//MuonRec.SetOption("MUON","SIMPLEFIT"); // new scheme simple fitting
//MuonRec.SetOption("MUON","COG"); // new scheme basic center-of-gravity only
// Use the following to disconnect the status map creation (only for debug!)
// as this speeds up startup a bit...
//MuonRec.SetOption("MUON","NOSTATUSMAP");
// Use the following to write to disk the digits (from raw data)
//MuonRec.SetOption("MUON","SAVEDIGITS");
MuonRec.Run();
.q
EOF

if [ ! -e MUON.Digits.root ]; then
    echo "Moving Digits files back ..."
    mv MUON.Digits/MUON.Digits.root .
fi 

echo "Running Trigger efficiency  ..."
aliroot -b >& testTriggerResults.out << EOF
.L $ALICE_ROOT/MUON/MUONTriggerEfficiency.C+
MUONTriggerEfficiency("galice.root",1);
.q
EOF

echo "Running efficiency  ..."

aliroot -b >& testResults.out << EOF
.L $ALICE_ROOT/MUON/MUONefficiency.C+
// no argument assumes Upsilon but MUONefficiency(443) works on Jpsi
MUONefficiency();
.q
EOF

echo "Running check ..."

aliroot -b >& testCheck.out << EOF
gSystem->Load("libMUONevaluation");
.L $ALICE_ROOT/MUON/MUONCheck.C+
MUONCheck(0, 9, "galice.root","AliESDs.root"); 
.q
EOF

echo "Running dumps for selected event (5) ..."

aliroot -b >& testDump.out << EOF
AliMUONData data("galice.root");
data.DumpKine(5);       > dump.kine
data.DumpHits(5);       > dump.hits
data.DumpDigits(5);     > dump.digits
data.DumpSDigits(5);    > dump.sdigits
data.DumpRecPoints(5);  > dump.recpoints
data.DumpRecTrigger(5); > dump.rectrigger
.q
EOF

echo "Finished"  
echo "... see results in test_out"

cd $CURDIR
