Int_t TestSPD(Int_t n = 5){
   Int_t rc=0;

   AliLoader::SetDebug(kTRUE);//set it to kTRUE for debug print-out
   gAlice->SetDebug(100);
   /**********************************************/
   /************ G E N E R A T I O N *************/
   /**********************************************/

   gROOT->LoadMacro("$(ALICE_ROOT)/macros/grun.C");
   grun(n,"ConfigSPD02.C");

   /**********************************************/
   /******************* I T S ********************/
   /**********************************************/
   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSHits2SDigits.C");
   AliITSHits2SDigits();
   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSSDigits2Digits.C");
   AliITSSDigits2Digits();
   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSDigits2RecPoints.C");
   AliITSDigits2RecPoints();

}
