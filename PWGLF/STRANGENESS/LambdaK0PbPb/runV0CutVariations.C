void runV0CutVariations(Double_t cmin, Double_t cmax, Bool_t isMC=kFALSE, Bool_t selectNonInjected=kFALSE) {
   TProof::Open("");

   TChain *ch=new TChain("chain");
   ch->SetProof();
   
   gProof->Load("AliV0CutVariations.C+");
   AliV0CutVariations *selector=new AliV0CutVariations();
   selector->SetCentrality(cmin,cmax);
   selector->SetMC(isMC);
   selector->SetSelectNonInjected(selectNonInjected);

   if (isMC) {
    ch->Add("LHC11a10b_plus/Merged.root/PWGLFExtractPerformanceV0_PP_MC/fTree");
   } else {
     ch->Add("LHC10h_pass2/Merged.root/PWGLFExtractV0_PP/fTree");
   }

   ch->Process(selector);

}
