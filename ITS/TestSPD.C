Int_t TestSPD(Int_t n = 10,Char_t SlowOrFast='s'){
   
   AliLoader::SetDebug(kTRUE);//set it to kTRUE for debug print-out
   gAlice->SetDebug(100);
   /**********************************************/
   /************ G E N E R A T I O N *************/
   /**********************************************/
      
   gROOT->LoadMacro("$(ALICE_ROOT)/macros/grun.C");
   grun(n,"$(ALICE_ROOT)/ITS/ConfigSPD02.C");

  
   /**********************************************/
   /******************* I T S ********************/
   /**********************************************/
   //The following part is just core of AliITStestV2.C
   if (SlowOrFast=='f') {
      cerr<<"Fast AliITSRecPoint(s) !\n";
      gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSHits2FastRecPoints.C");
      AliITSHits2FastRecPoints();
   } else {
      cerr<<"Slow AliITSRecPoint(s) !\n";
      gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSHits2SDigits.C");
      AliITSHits2SDigits();
      gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSSDigits2Digits.C");
      AliITSSDigits2Digits();
      gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSDigits2RecPoints.C");
      AliITSDigits2RecPoints();
   }
        // The following may not work properly.
//   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSFindClustersV2.C");
//   if (rc=AliITSFindClustersV2(SlowOrFast)) return rc;

//   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSFindTracksV2.C");
//   if (rc=AliITSFindTracksV2()) return rc;

   ::Info("NewIO test","Everything seems to be OK");   
   ::Info("NewIO test","You can try now display.C");
}
