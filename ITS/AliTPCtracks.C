Int_t AliTPCtracks() {
   Int_t rc=0;

//Test TPC simulation

   gROOT->LoadMacro("$(ALICE_ROOT)/macros/grun.C");
   grun();

   Int_t ver=gAlice->GetDetector("TPC")->IsVersion();
   delete gAlice; gAlice=0;

   if ((ver!=1)&&(ver!=2)) {
      cerr<<"Invalid TPC version: "<<ver<<" ! (must be 1 or 2)\n";
      return 12345;
   }

   if (ver==2) {
     gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCHits2Digits.C");
     if (rc=AliTPCHits2Digits()) return rc;

   }

//Test TPC reconstruction
   gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCFindClusters.C");
   if (rc=AliTPCFindClusters()) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCFindTracks.C");
   if (rc=AliTPCFindTracks()) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/TPCtracks.C");
   if (rc=TPCtracks()) return rc;

   return rc;
}
