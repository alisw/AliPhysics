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

CDBDIRECTORY="$ALICE_ROOT/MUON/CDB/Default";
CDB="local://$CDBDIRECTORY";

if [ ! -d $CDBDIRECTORY"/MUON" ]; then

echo "Generating Condition Database in directory $CDBDIRECTORY. This may take a while, so please be patient..."

aliroot -b >& testGenerateCalibrations.out << EOF
.L $ALICE_ROOT/MUON/MUONCDB.C++
gRandom->SetSeed($SEED);
generateCalibrations("$CDB",true);
.q
EOF

else

echo "Condition Database found in directory $CDBDIRECTORY. Will use it if needed."

fi

echo "Running simulation  ..."

aliroot -b  >& testSim.out << EOF 
gRandom->SetSeed($SEED);
AliCDBManager::Instance()->SetDefaultStorage("$CDB");
AliSimulation MuonSim("$ALICE_ROOT/MUON/Config.C");
MuonSim.SetWriteRawData("MUON");
MuonSim.Run($NEVENTS)
.q
EOF

echo "Running reconstruction  ..."

aliroot -b >& testReco.out << EOF
gRandom->SetSeed($SEED);
AliCDBManager::Instance()->SetDefaultStorage("$CDB");
AliReconstruction MuonRec("galice.root")
MuonRec.SetInput("$FULLPATH/");
MuonRec.SetRunTracking("")
MuonRec.SetRunVertexFinder(kFALSE)
MuonRec.SetRunLocalReconstruction("MUON")
MuonRec.SetFillESD("MUON")
MuonRec.Run();
.q
EOF

echo "Running Trigger efficiency  ..."
aliroot -b >& testTriggerResults.out << EOF
.L $ALICE_ROOT/MUON/MUONTriggerEfficiency.C++
MUONTriggerEfficiency();
.q
EOF

echo "Running efficiency  ..."

aliroot -b >& testResults.out << EOF
.L $ALICE_ROOT/MUON/MUONefficiency.C++
// no argument assumes Upsilon but MUONefficiency(443) works on Jpsi
MUONefficiency();
.q
EOF

echo "Running dumps ..."

if [ "$NEVENTS" -le 20 ]; then

aliroot -b > /dev/null << EOF
.L $ALICE_ROOT/MUON/MUONCheck.C++
MUONdigits(); > check.digits
MUONrecpoints(); > check.recpoints
MUONrectracks(); > check.rectracks
MUONrectrigger(); > check.rectrigger
EOF

fi

echo "Finished"  
echo "... see results in test_out"

cd $CURDIR
