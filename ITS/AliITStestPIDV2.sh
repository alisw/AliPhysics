#!/bin/sh
# 12 May 2003,Dubna
# This script processes the following steps for n events (AliRoot v3-09-08,
# root v3.05/04)
# (use the parameter n):
# - the TPC and ITS slow simulation (hits+digits+clusters),
# - the TPC+ITS tracking,
# - the TPC and ITS PID
# - the AliITStracksPidV2.root is created with the following information:
#    fGlab - track number
#    fPcode - particle code after the TPC PID
#    fMom   - particle momentum      (from the AliTPCtracks.root)
#    fLam   - particle lambder angle (from the AliTPCtracks.root)
#    fPhi   - particle phi angle     (from the AliTPCtracks.root)
#    fSignal- TPC signal (mips)
#    fWpi   - the PID weight for the identified pion
#    fWk    - the PID weight for the identified kaon
#    fWp    - the PID weight for the identified proton
#       (the weghts are tuned for the HIJING events)
# and the AliITSScanPIDV2.C is the example of macros to read the PID
# information.
# See ITSPIDplot.ps after run of this script.

if [ $# = 0   ]; then nev=1; else nev=$1; fi
if [ $nev = 0 ]; then nev=1; fi
#
# delete eventual old files from the last run
echo "Start simulation for " $nev " event(s)"
$ALICE_ROOT/ITS/AliITSDeleteOldFiles.sh
#
# run the hit generation
aliroot -q -b "$ALICE_ROOT/macros/grun.C($nev)"  
# digitize TPC
aliroot -q -b "$ALICE_ROOT/TPC/AliTPCHits2Digits.C($nev)" 
# TPC tracking
ln -s galice.root digits.root
aliroot -q -b "$ALICE_ROOT/TPC/AliTPCFindClustersMI.C($nev)" 
aliroot -b <<EOI
.L $ALICE_ROOT/TPC/AliTPCFindTracksMI.C
AliTPCFindTracks($nev);
.q
EOI
#
# digitize ITS
aliroot -q -b "$ALICE_ROOT/ITS/AliITSHits2SDigits.C"
aliroot -q -b "$ALICE_ROOT/ITS/AliITSSDigits2Digits.C"
# ITS-TPC tracking
aliroot -q -b "$ALICE_ROOT/ITS/AliITSFindClustersV2.C('s',$nev)"
aliroot -b <<EOI
TFile *inkin = TFile::Open("galice.root");
if (!inkin->IsOpen()) cerr<<"Can't open galice.root !\n";               
if (gAlice) {delete gAlice; gAlice=0;}
                                  
gAlice = (AliRun*)inkin->Get("gAlice");
cout<<"AliRun object found on file "<<gAlice<<endl;
cout<<"!!!! field ="<<gAlice->Field()->SolenoidField()<<endl;
AliKalmanTrack::SetConvConst(1000/0.299792458/gAlice->Field()->SolenoidField());
inkin->Close();
.x $ALICE_ROOT/ITS/AliITSFindTracksV2.C($nev);
.q
EOI
#
# Do the PID procedure for the ITS. The output file AliITStrackV2Pid.root is
# created with the information for each track:
aliroot -q -b "$ALICE_ROOT/ITS/AliITSSavePIDV2.C(0,$[$nev-1])"
aliroot -q -b "$ALICE_ROOT/ITS/AliITSScanPIDV2.C(0,$[$nev-1])"



