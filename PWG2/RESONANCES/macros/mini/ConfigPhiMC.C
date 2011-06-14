//
// *** Configuration script for phi->KK analysis with 2010 runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
//
Bool_t ConfigPhiMC(AliRsnMiniAnalysisTask *task, const char *suffix)
{
   // 
   // -- Define track cuts -------------------------------------------------------------------------
   //
   
   /*** EMPTY FOR TRUE PAIRS COMPUTATION ***/
   
   //
   // -- Values ------------------------------------------------------------------------------------
   //
   
   // invariant mass
   Int_t imID = task->CreateValue(AliRsnMiniValue::kInvMass, kTRUE);
   
   // transverse momentum
   Int_t ptID = task->CreateValue(AliRsnMiniValue::kPt, kTRUE);
   
   //
   // -- Pair cuts ---------------------------------------------------------------------------------
   //
   
   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidityMC", AliRsnCutMiniPair::kRapidityRangeMC);
   cutY->SetRangeD(-0.5, 0.5);
   
   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCutsMC", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   
   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput(Form("phi_TrueMC_%s", suffix), "HIST", "MOTHER");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kKaon);
   out->SetDaughter(1, AliRsnDaughter::kKaon);
   out->SetMotherPDG(333);
   out->SetMotherMass(1.019455);
   // pair cuts
   out->SetPairCuts(cutsPair);
   // binnings
   out->AddAxis(imID, 500, 0.9,  1.4);
   out->AddAxis(ptID, 200, 0.0, 10.0);
   
   return kTRUE;
}
