/****************************************************************************
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 ****************************************************************************/

Int_t AliTPCtest(Int_t n = 5) {
   Int_t rc=0;

//Test TPC simulation
   gROOT->LoadMacro("$(ALICE_ROOT)/macros/grun.C");
   grun(n);

   
AliKalmanTrack::SetConvConst(1000/0.299792458/gAlice->Field()->SolenoidField());

   Int_t ver=gAlice->GetDetector("TPC")->IsVersion();
   
   AliRunLoader *rl = gAlice->GetRunLoader();
   if (rl == 0x0) {
      cerr<<"Can not get run loader from gAlice"<<endl;
      return 1;
   }

   delete rl;//close the session left after generation (grun.C)
   delete gAlice; gAlice=0x0;
   
   cout<<" \n\n\nClean -> Proceeding witg digitization \n\n\n";
   if ((ver!=1)&&(ver!=2)) {
      cerr<<"Invalid TPC version: "<<ver<<" ! (must be 1 or 2)\n";
      return 12345;
   }

   if (ver==2) {
     gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCHits2Digits.C");
     if (rc=AliTPCHits2Digits(n)) return rc;

     //gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCDisplayDigits.C");
     //if (rc=AliTPCDisplayDigits(1,1)) return rc;
   }


//Test TPC reconstruction
   gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCFindClusters.C");
   if (rc=AliTPCFindClusters(n)) return rc;

   //gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCDisplayClusters.C");
   //if (rc=AliTPCDisplayClusters()) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCFindTracks.C");
   if (rc=AliTPCFindTracks(n)) return rc;

   gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS");
   gROOT->ProcessLine(".L $(ALICE_ROOT)/TPC/AliTPCComparison.C++");
   if (rc=AliTPCComparison()) return rc;

   return rc;
}

