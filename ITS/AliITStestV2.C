/****************************************************************************
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 ****************************************************************************/

Int_t AliITStestV2(Int_t nev=5, Char_t SlowOrFast='s') {
   Int_t rc=0;

   if (gAlice) {
      delete gAlice->GetRunLoader();
      delete gAlice; 
      gAlice=0;
   }

   rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0) {
      cerr<<"Can not open session"<<endl;
      return 1;
   }

   if (rl->LoadgAlice()) {
      cerr<<"Error occured while loading AliRun"<<endl;
      return 1;
   }
   AliKalmanTrack::
   SetConvConst(1000/0.299792458/rl->GetAliRun()->Field()->SolenoidField());

   delete rl;
  
   if (SlowOrFast=='f') {
      cerr<<"Fast AliITSRecPoint(s) !\n";
   } else {
      cerr<<"Slow AliITSRecPoint(s) !\n";
      gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSHits2SDigits.C");
      AliITSHits2SDigits();
      gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSSDigits2Digits.C");
      AliITSSDigits2Digits();
   }
   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSFindClustersV2.C");
   if (rc=AliITSFindClustersV2(SlowOrFast)) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSFindTracksV2.C");
   if (rc=AliITSFindTracksV2()) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSComparisonV2.C");
   if (rc=AliITSComparisonV2()) return rc;
/*
   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliV0FindVertices.C");
   if (rc=AliV0FindVertices()) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliV0Comparison.C");
   if (rc=AliV0Comparison()) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliCascadeFindVertices.C");
   if (rc=AliCascadeFindVertices()) return rc;

   gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS");
   gROOT->ProcessLine(".L $(ALICE_ROOT)/ITS/AliCascadeComparison.C+");
   if (rc=AliCascadeComparison()) return rc;
*/
   return rc;
}
