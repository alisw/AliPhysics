void runV0CutVariations(Bool_t isMC=kFALSE) {
   TProof::Open("");

   TChain *ch=new TChain("chain");
   ch->SetProof();
   
   gProof->Load("AliV0CutVariations.C+");
   AliV0CutVariations *selector=new AliV0CutVariations();
   selector->SetMC(isMC);

   if (isMC) {
     ch->Add("DavidsV0MC_offline.root/PWG2CheckPerformanceLambda_PP_MC/fTree");
   } else {
     ch->AddFile("DavidsV0_off.root/PWG2CheckLambda_off_PP/fTree");
   }

   ch->Process(selector);

}
