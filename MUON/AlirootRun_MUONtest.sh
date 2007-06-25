#!/bin/sh
# $Id$
# with galice.root, galice_sim.root 

if [ $# -ne 1 ]; then
  NEVENTS=100
  echo "Number of events not specified. Using $NEVENTS"
else
  NEVENTS=$1
fi

CURDIR=`pwd`
OUTDIR=test_out.$NEVENTS

rm -fr $OUTDIR
mkdir $OUTDIR
cp $ALICE_ROOT/MUON/.rootrc $ALICE_ROOT/MUON/rootlogon.C $OUTDIR
cd $OUTDIR

DUMPEVENT=5
RUN=0
FULLPATH="$CURDIR/$OUTDIR"
SEED=1234567
SIMDIR="generated"

echo "Running simulation  ..."

aliroot -b >& testSim.out << EOF 
// Uncoment following lines to run simulation with local residual mis-alignment
// (generated via MUONGenerateGeometryData.C macro)
// AliCDBManager* man = AliCDBManager::Instance();
// man->SetDefaultStorage("local://$ALICE_ROOT");
// man->SetSpecificStorage("MUON/Align/Data","local://$ALICE_ROOT/MUON/ResMisAlignCDB");
gRandom->SetSeed($SEED);
AliCDBManager::Instance()->SetRun($RUN);
AliSimulation MuonSim("$ALICE_ROOT/MUON/Config.C");
MuonSim.SetMakeTrigger("MUON");
MuonSim.SetWriteRawData("MUON","raw.root",kTRUE);
MuonSim.Run($NEVENTS);
.q
EOF

echo "Moving generated files to $SIMDIR"
mkdir $SIMDIR
mv MUON*.root Kinematics*.root galice.root TrackRefs*.root $SIMDIR

echo "Running reconstruction  ..."

aliroot -b >& testReco.out << EOF
AliCDBManager::Instance()->SetRun($RUN);
gRandom->SetSeed($SEED);
AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k5kG);
AliTracker::SetFieldMap(field, kFALSE);
AliReconstruction* MuonRec = new AliReconstruction("galice.root");
MuonRec.SetInput("$FULLPATH/raw.root");
MuonRec.SetRunVertexFinder(kFALSE);
MuonRec.SetRunLocalReconstruction("MUON");
MuonRec.SetRunTracking("MUON");
MuonRec.SetFillESD("");
MuonRec.SetLoadAlignData("MUON");
MuonRec.SetNumberOfEventsPerFile(1000);
MuonRec.SetOption("MUON","SAVEDIGITS");
// Change the line above with the line below to run without local reconstruction
// MuonRec.SetOption("MUON","SAVEDIGITS NOLOCALRECONSTRUCTION");
MuonRec.Run();
delete MuonRec;
.q
EOF

echo "Running Trigger efficiency  ..."
aliroot -b >& testTriggerResults.out << EOF
.L $ALICE_ROOT/MUON/MUONTriggerEfficiency.C+
MUONTriggerEfficiency("$SIMDIR/galice.root", "galice.root", 1);
.q
EOF

echo "Running efficiency  ..."

aliroot -b >& testResults.out << EOF
.L $ALICE_ROOT/MUON/MUONefficiency.C+
// no argument assumes Upsilon but MUONefficiency(443) works on Jpsi
MUONefficiency("$SIMDIR/galice.root");
.q
EOF

echo "Running check ..."
aliroot -b >& testCheck.out << EOF
gSystem->Load("libMUONevaluation");
.L $ALICE_ROOT/MUON/MUONCheck.C+
MUONCheck(0, $NEVENTS-1, "generated/galice.root", "galice.root", "AliESDs.root"); 
.q
EOF


echo "Running dumps for selected event ($DUMPEVENT) ..."

aliroot -b  << EOF
AliMUONMCDataInterface mcdSim("$SIMDIR/galice.root");
mcdSim.DumpKine($DUMPEVENT);       > dump.kine
mcdSim.DumpHits($DUMPEVENT);       > dump.hits
mcdSim.DumpTrackRefs($DUMPEVENT);  > dump.trackrefs
mcdSim.DumpDigits($DUMPEVENT,true);  > dump.simdigits
mcdSim.DumpSDigits($DUMPEVENT,true); > dump.sdigits

AliMUONDataInterface dRec("galice.root");
dRec.DumpRecPoints($DUMPEVENT);  > dump.recpoints
dRec.DumpTracks($DUMPEVENT);     > dump.tracks
dRec.DumpTriggerTracks($DUMPEVENT); > dump.triggertracks
dRec.DumpTrigger($DUMPEVENT);  > dump.trigger
dRec.DumpDigits($DUMPEVENT,true); > dump.recdigits
.q
EOF


echo "Finished"  
echo "... see results in test_out"

cd $CURDIR
