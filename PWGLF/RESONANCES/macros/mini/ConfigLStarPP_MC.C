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
Bool_t ConfigLStarPP_MC
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
   
   TString mode = "SPARSE";
//   if (!isPP) mode = "SPARSE";
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput(Form("LStar_TrueMC1%s", suffix), mode.Data(), "MOTHER");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kProton);
   out->SetDaughter(1, AliRsnDaughter::kKaon);
      out->SetCharge(0, '+');
      out->SetCharge(1, '-');
   out->SetMotherPDG(3124);
   out->SetMotherMass(1.520);
   // pair cuts
   out->SetPairCuts(cutsPair);
   // binnings
   out->AddAxis(imID,  80, 1.0, 2.0);
//S.K.   out->AddAxis(ptID, 100, 0.0, 10.0);
   out->AddAxis(ptID, 100, 0.0, 10.0);
//S.K.   if (!isPP) out->AddAxis(centID, 100, 0.0, 100.0);
   out->AddAxis(centID, 10, 0.0, 100.0);
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput(Form("LStar_TrueMC2%s", suffix), mode.Data(), "MOTHER");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kProton);
   out->SetDaughter(1, AliRsnDaughter::kKaon);
   out->SetCharge(0, '-');
   out->SetCharge(1, '+');
   out->SetMotherPDG(-3124);
   out->SetMotherMass(1.520);
   // pair cuts
   out->SetPairCuts(cutsPair);
   // binnings
   out->AddAxis(imID,  100, 1.0, 2.0);
//S.K.   out->AddAxis(ptID, 100, 0.0, 10.0);
   out->AddAxis(ptID, 100, 0.0, 10.0);
//S.K.   if (!isPP) out->AddAxis(centID, 100, 0.0, 100.0);
   out->AddAxis(centID, 10, 0.0, 100.0);
     
    // create output
   AliRsnMiniOutput *outm = task->CreateOutput(Form("Ls_Mother"), mode.Data(), "MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kProton);
    outm->SetDaughter(1, AliRsnDaughter::kKaon);
    outm->SetMotherPDG(3124);
    outm->SetMotherMass(1.520);
    // pair cuts
    outm->SetPairCuts(cutsPair);
    // binnings
    outm->AddAxis(imID, 100, 1.0, 2.0);
    outm->AddAxis(ptID, 100, 0.0, 10.0);
    if (!isPP){
      outm->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outm->AddAxis(centID, 400, 0.0, 400.0);
    }
  
 AliRsnMiniOutput *outm = task->CreateOutput(Form("Ls_AntiMother"), mode.Data(), "MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kProton);
    outm->SetDaughter(1, AliRsnDaughter::kKaon);
    outm->SetMotherPDG(-3124);
    outm->SetMotherMass(1.520);
    // pair cuts
    outm->SetPairCuts(cutsPair);
    // binnings
    outm->AddAxis(imID, 100, 1.0, 2.0);
    outm->AddAxis(ptID, 100, 0.0, 10.0);
    if (!isPP){
      outm->AddAxis(centID, 100, 0.0, 100.0);
    }   else    { 
      outm->AddAxis(centID, 400, 0.0, 400.0);
    }
   return kTRUE;
}
