Int_t AliITStestV2(Char_t SlowOrFast='s') {
   Int_t rc=0;

   if (gAlice) {delete gAlice; gAlice=0;}
   TFile *in=TFile::Open("galice.root");
   if (!in->IsOpen()) {
      cerr<<"Can't open galice.root !\n"; 
      return 1;
   }
   if (!(gAlice=(AliRun*)in->Get("gAlice"))) {
      cerr<<"Can't find gAlice !\n";
      return 2;
   }
   AliKalmanTrack::SetConvConst(100/0.299792458/0.2/gAlice->Field()->Factor());
   delete gAlice; gAlice=0;
   in->Close();

   if (SlowOrFast=='f') {
      cerr<<"Fast AliITSRecPoint(s) !\n";
      gROOT->LoadMacro("$(ALICE_ROOT)/ITS/ITSHitsToFastPoints.C");
      ITSHitsToFastPoints();
   } else {
      cerr<<"Slow AliITSRecPoint(s) !\n";
      gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSHits2Digits.C");
      AliITSHits2Digits();
      gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSFindClusters.C");
      AliITSFindClusters();
   }

   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSFindClustersV2.C");
   if (rc=AliITSFindClustersV2()) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSFindTracksV2.C");
   if (rc=AliITSFindTracksV2()) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSComparisonV2.C");
   if (rc=AliITSComparisonV2()) return rc;

   return rc;
}
