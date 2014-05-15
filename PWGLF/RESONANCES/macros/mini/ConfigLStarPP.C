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
Bool_t ConfigLStarPP
(  
   AliRsnMiniAnalysisTask *task, 
   Bool_t                  isMC, 
   const char             *suffix,
   AliRsnCutSet           *cutsPair
)
{
   // manage suffix
   if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
   
   // 
   // -- Define track cuts -------------------------------------------------------------------------
   //
   
   // integrated pion cut
   AliRsnCutDaughterLStar2010 *cutK = new AliRsnCutDaughterLStar2010("cutKaonForLStar", AliPID::kKaon);
   // cut set
   AliRsnCutSet *cutSetK = new AliRsnCutSet("setPionForLStar", AliRsnTarget::kDaughter);
   cutSetK->AddCut(cutK);
   cutSetK->SetCutScheme(cutK->GetName());
   // add to task
   Int_t icutK = task->AddTrackCuts(cutSetK);
   
   // integrated kaon cut
   AliRsnCutDaughterLStar2010 *cutP = new AliRsnCutDaughterLStar2010("cutProtonForLStar", AliPID::kProton);
   // cut set
   AliRsnCutSet *cutSetP = new AliRsnCutSet("setKaonForLStar", AliRsnTarget::kDaughter);
   cutSetP->AddCut(cutP);
   cutSetP->SetCutScheme(cutP->GetName());
   // add to task
   Int_t icutP = task->AddTrackCuts(cutSetP);
   
   //
   // -- Values ------------------------------------------------------------------------------------
   //
   
   /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
   /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
   /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
   /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
   
   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //
   
   // use an array for more compact writing, which are different on mixing and charges
   // [0] = unlike
   // [1] = mixing
   // [2] = like ++
   // [3] = like --
   Bool_t  use     [10] = { 1       ,  1       ,  1       ,  1       ,  1      ,  1      ,  isMC   ,   isMC   ,  isMC   ,   isMC   };
   Bool_t  useIM   [10] = { 1       ,  1       ,  1       ,  1       ,  1      ,  1      ,  1      ,   1      ,  0      ,   0      };
   TString name    [10] = {"Unlike1", "Unlike2", "Mixing1", "Mixing2", "LikePP", "LikeMM", "Trues1",  "Trues2", "Res1"  ,  "Res2"  };
   TString comp    [10] = {"PAIR"   , "PAIR"   , "MIX"    , "MIX"    , "PAIR"  , "PAIR"  , "TRUE"  ,  "TRUE"  , "TRUE"  ,  "TRUE"  };
   //TString output  [10] = {"HIST" , "HIST"   , "HIST"   , "HIST"   , "HIST"  , "HIST"  , "HIST"  ,  "HIST"  , "HIST"  ,  "HIST"  };
   TString output  [10] = {"SPARSE" , "SPARSE" , "SPARSE" , "SPARSE" , "SPARSE", "SPARSE", "SPARSE",  "SPARSE", "SPARSE",  "SPARSE"};
   Int_t   pdgCode [10] = { 3124    ,  -3124    ,  3124    , -3124    ,  3124   ,  3124   ,  3124   ,  -3124   ,  3124   ,  -3124  };
   Char_t  charge1 [10] = {'+'      , '-'      , '+'      , '-'      , '+'     , '-'     , '+'     ,  '-'     , '+'     ,  '-'     };
   Char_t  charge2 [10] = {'-'      , '+'      , '-'      , '+'      , '+'     , '-'     , '-'     ,  '+'     , '-'     ,  '+'     };
   Int_t   cutID1  [10] = { icutP   ,  icutP   ,  icutP   ,  icutP   ,  icutP  ,  icutP  ,  icutP  ,   icutP  ,  icutP  ,   icutP  };
   Int_t   cutID2  [10] = { icutK   ,  icutK   ,  icutK   ,  icutK   ,  icutK  ,   icutK ,   icutK ,   icutK  ,  icutK  ,   icutK  };
   for (Int_t i = 0; i < 10; i++) {
      if (!use[i]) continue;
      // create output
      AliRsnMiniOutput *out = task->CreateOutput(Form("LStar_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      out->SetCutID(0, cutID1[i]);
      out->SetCutID(1, cutID2[i]);
      out->SetDaughter(0, AliRsnDaughter::kProton);
      out->SetDaughter(1, AliRsnDaughter::kKaon);
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetMotherPDG(pdgCode[i]);
      out->SetMotherMass(1.520);
      // pair cuts
      out->SetPairCuts(cutsPair);
      // axis X: invmass (or resolution)
      if (useIM[i]) 
         out->AddAxis(imID, 80, 1.4, 1.8);
      else
         out->AddAxis(resID, 200, -0.02, 0.02);
      // axis Y: transverse momentum
//S.K. out->AddAxis(ptID, 100, 0.0, 10.0);
      out->AddAxis(ptID, 100, 0.0, 10.0);
      //S.K. axis Z: centrality 
      out->AddAxis(centID, 10, 0.0, 100.0);
   }
   
   return kTRUE;
}
