Int_t AliITStestBari() {
   Int_t rc=0;

//Test ITS simulation
   gROOT->LoadMacro("$(ALICE_ROOT)/macros/grun.C");
   grun();

   Int_t ver=gAlice->GetDetector("ITS")->IsVersion();
   delete gAlice; gAlice=0;

   if (ver!=5) {
      cerr<<"Invalid ITS version: "<<ver<<" ! (must be 5 for the moment)\n";
      return 12345;
   }

   if (ver==5) {
     gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSHits2DigitsBari.C");
     if (rc=AliITSHits2Digits()) return rc;

   }

   printf("start reconstruction\n");

//Test ITS reconstruction
   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSFindClustersBari.C");
   if (rc=AliITSFindClusters()) return rc;

   //gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSgraphycs.C");
   //if (rc=AliITSgraphycs()) return rc;

   return rc;
}
