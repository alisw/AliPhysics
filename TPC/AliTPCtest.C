Int_t AliTPCtest() {
   Int_t rc=0;

//Test TPC simulation
   gROOT->LoadMacro("$(ALICE_ROOT)/macros/grun.C");
   grun();

   Int_t ver=gAlice->GetDetector("TPC")->IsVersion();
   delete gAlice; gAlice=0;
   if (ver==2) {
     gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCHits2Digits.C");
     if (rc=AliTPCHits2Digits()) return rc;

     gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCDisplayDigits.C");
     if (rc=AliTPCDisplayDigits(1,1)) return rc;
   }

//Test TPC reconstruction
   gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCFindClusters.C");
   if (rc=AliTPCFindClusters()) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCDisplayClusters.C");
   if (rc=AliTPCDisplayClusters()) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCFindTracks.C");
   if (rc=AliTPCFindTracks()) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCComparison.C");
   if (rc=AliTPCComparison()) return rc;

   return rc;
}
