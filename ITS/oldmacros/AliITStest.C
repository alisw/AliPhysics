Int_t AliITStest() {
   Int_t rc=0;

//Test ITS simulation
   gROOT->LoadMacro("$(ALICE_ROOT)/macros/grun.C");
   grun();

   AliITSgeom *gm = ((AliITS*)(gAlice->GetDetector("ITS")))->GetITSgeom();
   delete gAlice; gAlice=0;

   if (!gm) {
       cerr << "This version of the ITS geometry does not have a"
	    << " AliITSgeom defined" << endl;
      return 12345;
   }

   if (gm) {
     gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSHits2DigitsDefault.C");
     if (rc=AliITSHits2DigitsDefault()) return rc;
   }

   cout << "start reconstruction" << endl;

   //Test ITS reconstruction
   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSFindClusters.C");

   delete gAlice; gAlice=0;
   
   if (rc=AliITSFindClusters()) return rc;

   //gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSgraphycs.C");
   //if (rc=AliITSgraphycs()) return rc;

   return rc;
}
