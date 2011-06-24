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
Bool_t ConfigSigmaStarMC
(
   AliRsnMiniAnalysisTask *task, 
   Bool_t                  isPP, 
   const char             *suffix,
   AliRsnCutSet           *cutsPair
)
{
   // manage suffix
   if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
   
   // 
   // -- Define track cuts -------------------------------------------------------------------------
   //
   
   /*** EMPTY FOR TRUE PAIRS COMPUTATION ***/
   
   //
   // -- Values ------------------------------------------------------------------------------------
   //
   
   /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
   /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
   /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);

   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //
   
   TString mode = "HIST";
   if (!isPP) mode = "SPARSE";
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarP_TrueMC%s", suffix), mode.Data(), "MOTHER");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kLambda);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetMotherPDG(3224);
   out->SetMotherMass(1382.3);
   // pair cuts
   out->SetPairCuts(cutsPair);
   // binnings
   out->AddAxis(imID, 800, 1.2, 2.0);
   out->AddAxis(ptID, 100, 0.0, 10.0);
   if (!isPP) out->AddAxis(centID, 100, 0.0, 100.0);
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarM_TrueMC%s", suffix), mode.Data(), "MOTHER");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kLambda);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetMotherPDG(3114);
   out->SetMotherMass(1387.4);
   // pair cuts
   out->SetPairCuts(cutsPair);
   // binnings
   out->AddAxis(imID, 800, 1.2, 2.0);
   out->AddAxis(ptID, 100, 0.0, 10.0);
   if (!isPP) out->AddAxis(centID, 100, 0.0, 100.0);
   
   return kTRUE;
}
