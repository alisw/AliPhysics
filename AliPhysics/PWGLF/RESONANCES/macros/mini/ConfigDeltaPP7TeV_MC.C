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
Bool_t ConfigDeltaPP7TeV_MC
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
   /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt     , kFALSE);
   /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult   , kFALSE);
      
   /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta , kFALSE);
   /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY   , kFALSE); 
   
   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //
   
   TString mode = "SPARSE";
   //if (!isPP) mode = "SPARSE";
   
   
   //DELTA++
   // create output
   AliRsnMiniOutput *outpp = task->CreateOutput(Form("deltapp_TrueMCGen%s", suffix), mode.Data(), "MOTHER");
   outpp->SetDaughter(0, AliRsnDaughter::kProton);
   outpp->SetDaughter(1, AliRsnDaughter::kPion);
   outpp->SetMotherPDG(2224);
   outpp->SetMotherMass(1.232);
   outpp->SetPairCuts(cutsPair);
   outpp->AddAxis( imID   , 100 , 1.0, 2.0  );
   outpp->AddAxis( ptID   , 100 , 0.0, 10.0 );
   outpp->AddAxis(centID  , 120 , 0.0, 120  );
   outpp->AddAxis( etaID  , 20  ,-1.0, 1.0  );
   outpp->AddAxis( yID    , 10  ,-0.5, 0.5  ); 
   

   //DELTA--
   AliRsnMiniOutput *outmm = task->CreateOutput(Form("deltamm_TrueMCGen%s", suffix), mode.Data(), "MOTHER");
   outmm->SetDaughter(0, AliRsnDaughter::kProton);
   outmm->SetDaughter(1, AliRsnDaughter::kPion);
   outmm->SetMotherPDG(-2224);
   outmm->SetMotherMass(1.232);
   outmm->SetPairCuts(cutsPair);
   outmm->AddAxis( imID   , 100 , 1.0, 2.0  );
   outmm->AddAxis( ptID   , 100 , 0.0, 10.0 );
   outmm->AddAxis( centID , 120 , 0.0, 120  );
   outmm->AddAxis( etaID  , 20  ,-1.0, 1.0  );
   outmm->AddAxis( yID    , 10  ,-0.5, 0.5  ); 
   
   
   //DELTA0
   AliRsnMiniOutput *out0pm = task->CreateOutput(Form("delta0pm_TrueMCGen%s", suffix), mode.Data(), "MOTHER");
   out0pm->SetDaughter(0, AliRsnDaughter::kProton);
   out0pm->SetDaughter(1, AliRsnDaughter::kPion);
   out0pm->SetMotherPDG(2114);
   out0pm->SetMotherMass(1.232);
   out0pm->SetPairCuts(cutsPair); 
   out0pm->AddAxis( imID   , 100 , 1.0, 2.0  );
   out0pm->AddAxis( ptID   , 100 , 0.0, 10.0 );
   out0pm->AddAxis( centID , 120 , 0.0, 120  );
   out0pm->AddAxis( etaID  , 20  ,-1.0, 1.0  );
   out0pm->AddAxis( yID    , 10  ,-0.5, 0.5  ); 
   
   //DELTA0BAR
   AliRsnMiniOutput *out0mp = task->CreateOutput(Form("delta0mp_TrueMCGen%s", suffix), mode.Data(), "MOTHER");
   out0mp->SetDaughter(0, AliRsnDaughter::kProton);
   out0mp->SetDaughter(1, AliRsnDaughter::kPion);
   out0mp->SetMotherPDG(-2114);
   out0mp->SetMotherMass(1.232);
   out0mp->SetPairCuts(cutsPair);
   out0mp->AddAxis( imID   , 100 , 1.0, 2.0  );
   out0mp->AddAxis( ptID   , 100 , 0.0, 10.0 );
   out0mp->AddAxis( centID , 120 , 0.0, 120  );
   out0mp->AddAxis( etaID  , 20  ,-1.0, 1.0  );
   out0mp->AddAxis( yID    , 10  ,-0.5, 0.5  ); 
   
   
   return kTRUE;
}
