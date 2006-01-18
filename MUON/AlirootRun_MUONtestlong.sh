#!/bin/sh
# $Id$

CURDIR=`pwd`
OUTDIR=testlong_out

rm -fr $OUTDIR
mkdir $OUTDIR
cp .rootrc $OUTDIR
cd $OUTDIR

echo "Running simulation  ..."

aliroot -b >& testSim.out << EOF  
AliSimulation MuonSim
MuonSim.SetConfigFile("$ALICE_ROOT/MUON/Config.C")
// Minimum number of events to have enough stat. for invariant mass fit
// 10000 is ok, 20000 is really fine
MuonSim.Run(1000) 
.q
EOF

echo "Running reconstruction  ..."

aliroot -b >& testReco.out << EOF 
TPluginManager* pluginManager = gROOT->GetPluginManager();
pluginManager->AddHandler("AliReconstructor", "MUON","AliMUONReconstructor", "MUON","AliMUONReconstructor()")
AliReconstruction MuonRec("galice.root") 
MuonRec.SetRunTracking("")
MuonRec.SetRunVertexFinder(kFALSE)
MuonRec.SetRunLocalReconstruction("MUON")
MuonRec.SetFillESD("MUON")
MuonRec.Run() 
.q
EOF

echo "Running Trigger efficiency  ..."

aliroot -b >& testTriggerResults.out << EOF
.includepath $ALICE_ROOT/STEER
.includepath $ALICE_ROOT/MUON
.L $ALICE_ROOT/MUON/MUONTriggerEfficiency.C++
MUONTriggerEfficiency();
.q
EOF

echo "Running efficiency  ..."

aliroot -b >& testEfficiency.out << EOF 
.includepath $ALICE_ROOT/STEER
.includepath $ALICE_ROOT/MUON
.L $ALICE_ROOT/MUON/MUONefficiency.C++
// no argument assumes Upsilon but MUONefficiency(443) handles Jpsi
MUONefficiency();
.q
EOF


aliroot -b >& testResults.out << EOF 
// no argument assumes Upsilon but MUONplotefficiency(443) handles Jpsi
.x $ALICE_ROOT/MUON/MUONplotefficiency.C
.q
EOF

more  testSim.out | grep 'RunSimulation: Execution time:'  > testTime.out
more  testSim.out | grep 'RunSDigitization: Execution time:'  >> testTime.out
more  testSim.out | grep 'RunDigitization: Execution time:'  >> testTime.out 

more  testReco.out | grep 'RunLocalReconstruction: Execution time for MUON'  >> testTime.out
more  testReco.out | grep 'Execution time for filling ESD ' >> testTime.out

rm gphysi.dat
rm *.root
rm testSim.out
rm testReco.out
rm *.eps

echo "Finished"  
echo "... see results in testlong_out"

cd $CURDIR
