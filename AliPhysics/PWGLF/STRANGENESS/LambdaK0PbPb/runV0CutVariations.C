void runV0CutVariations(
Double_t cmin, Double_t cmax, 
Bool_t isMC=kFALSE, Bool_t selectNonInjected=kFALSE,
TString fname) {
   TProof::Open("");

   TChain *ch=new TChain("chain");
   ch->SetProof();
   
   gProof->Load("AliV0CutVariations.C+");
   AliV0CutVariations *selector=new AliV0CutVariations();
   selector->SetCentrality(cmin,cmax);
   selector->SetMC(isMC);
   selector->SetSelectNonInjected(selectNonInjected);

   if (isMC) {
      fname += "/PWGLFExtractPerformanceV0_PP_MC/fTree";
   } else {
      fname += "/PWGLFExtractV0_PP/fTree";
   }

   cout<<"Running over "<<fname.Data()<<endl;

   ch->Add(fname.Data());  
   ch->Process(selector);
}
