void runV0CutVariations(Bool_t isMC=kTRUE, Bool_t SelectNonInjected=kFALSE) {
   TProof::Open("");

   TChain *ch=new TChain("chain");
   ch->SetProof();
   
   gProof->Load("AliV0CutVariations.C+");
   AliV0CutVariations *selector=new AliV0CutVariations();
   selector->SetMC(isMC);
   selector->SetSelectNonInjected(SelectNonInjected);

   if (isMC) {
   //ch->Add("LHC11a10a_bis/Merged.root/PWGLFExtractPerformanceV0_PP_MC/fTree");
    ch->Add("LHC11a10b_plus/Merged.root/PWGLFExtractPerformanceV0_PP_MC/fTree");
   } else {
     ch->Add("LHC10h_pass2/Merged.root/PWGLFExtractPerformanceV0_PP_MC/fTree");
   }

   ch->Process(selector);

}
