Int_t AliTPCtest() {
   Int_t rc=0;

//Test TPC simulation
   gROOT->LoadMacro("$(ALICE_ROOT)/macros/grun.C");
   grun();

   

AliKalmanTrack::SetConvConst(1000/0.299792458/gAlice->Field()->SolenoidField());

   Int_t ver=gAlice->GetDetector("TPC")->IsVersion();
   delete gAlice; gAlice=0;

   if ((ver!=1)&&(ver!=2)) {
      cerr<<"Invalid TPC version: "<<ver<<" ! (must be 1 or 2)\n";
      return 12345;
   }

   if (ver==2) {
     gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCHits2Digits.C");
     if (rc=AliTPCHits2Digits()) return rc;

     //     gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCDisplayDigits.C");
     //     if (rc=AliTPCDisplayDigits(1,1)) return rc;
   }


//Test TPC reconstruction
   gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCFindClusters.C");
   if (rc=AliTPCFindClusters()) return rc;

   //  gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCDisplayClusters.C");
   // if (rc=AliTPCDisplayClusters()) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCFindTracks.C");
   if (rc=AliTPCFindTracks()) return rc;

   gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS");
   gROOT->ProcessLine(".L $(ALICE_ROOT)/TPC/AliTPCComparison.C+");
   if (rc=AliTPCComparison()) return rc;

   return rc;
}

