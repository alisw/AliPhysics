#!/bin/sh
# $Id$

CURDIR=`pwd`
OUTDIR=testlong_out

rm -fr $OUTDIR
mkdir $OUTDIR
cp $ALICE_ROOT/MUON/.rootrc $ALICE_ROOT/MUON/rootlogon.C $OUTDIR
cd $OUTDIR

# Minimum number of events to have enough stat. for invariant mass fit
# 10000 is ok, 20000 is really fine
NEVENTS=10000
SEED=1234567


echo "Running simulation  ..."

aliroot -b >& testSim.out << EOF  
// Uncoment following lines to run simulation with local residual mis-alignment
// (generated via MUONGenerateGeometryData.C macro)
// AliCDBManager* man = AliCDBManager::Instance();
// man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
// man->SetSpecificStorage("MUON/Align/Data","local://$ALICE_ROOT/OCDB/MUON/ResMisAlignCDB");
AliSimulation MuonSim("$ALICE_ROOT/MUON/Config.C");
MuonSim.SetSeed($SEED);
MuonSim.SetMakeTrigger("MUON");
MuonSim.Run($NEVENTS); 
.q
EOF

echo "Running reconstruction  ..."

aliroot -b >& testReco.out << EOF 
man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
gRandom->SetSeed($SEED);
AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k5kG);
AliTracker::SetFieldMap(field, kFALSE);
AliReconstruction MuonRec("galice.root"); 
MuonRec.SetRunVertexFinder(kFALSE);
MuonRec.SetRunLocalReconstruction("MUON");
MuonRec.SetRunTracking("MUON");
MuonRec.SetFillESD("MUON");
MuonRec.SetLoadAlignData("MUON");
MuonRec.SetNumberOfEventsPerFile($NEVENTS);
MuonRec.Run(); 
.q
EOF

echo "Running Trigger efficiency  ..."

aliroot -b >& testTriggerResults.out << EOF
.L $ALICE_ROOT/MUON/MUONTriggerEfficiency.C+
MUONTriggerEfficiency("galice.root", "galice.root",0);
.q
EOF

echo "Running efficiency  ..."

aliroot -b >& testEfficiency.out << EOF 
.L $ALICE_ROOT/MUON/MUONefficiency.C+
// no argument assumes Upsilon but MUONefficiency(443) handles Jpsi
MUONefficiency("galice.root");
.q
EOF


aliroot -b >& testResults.out << EOF 
// no argument assumes Upsilon but MUONplotefficiency(443) handles Jpsi
.L $ALICE_ROOT/MUON/MUONplotefficiency.C+
MUONplotefficiency();
.q
EOF

more  testSim.out | grep 'RunSimulation: Execution time:'  > testTime.out
more  testSim.out | grep 'RunSDigitization: Execution time:'  >> testTime.out
more  testSim.out | grep 'RunDigitization: Execution time:'  >> testTime.out 

more  testReco.out | grep 'RunLocalReconstruction: Execution time for MUON'  >> testTime.out
more  testReco.out | grep 'Execution time for filling ESD ' >> testTime.out

#rm gphysi.dat
#rm *.root
#rm testSim.out
#rm testReco.out
#rm *.eps

echo "Finished"  
echo "... see results in testlong_out"

cd $CURDIR
