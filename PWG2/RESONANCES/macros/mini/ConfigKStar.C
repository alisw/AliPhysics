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
Bool_t ConfigKStar
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
   AliRsnCutDaughterKStar2010PP *cutPi = new AliRsnCutDaughterKStar2010PP("cutPionForKStar", AliPID::kPion);
   // cut set
   AliRsnCutSet *cutSetPi = new AliRsnCutSet("setPionForKStar", AliRsnTarget::kDaughter);
   cutSetPi->AddCut(cutPi);
   cutSetPi->SetCutScheme(cutPi->GetName());
   // add to task
   Int_t iCutPi = task->AddTrackCuts(cutSetPi);
   
   // integrated kaon cut
   AliRsnCutDaughterKStar2010PP *cutK = new AliRsnCutDaughterKStar2010PP("cutKaonForKStar", AliPID::kKaon);
   // cut set
   AliRsnCutSet *cutSetK = new AliRsnCutSet("setKaonForKStar", AliRsnTarget::kDaughter);
   cutSetK->AddCut(cutK);
   cutSetK->SetCutScheme(cutK->GetName());
   // add to task
   Int_t iCutK = task->AddTrackCuts(cutSetK);
   
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
   TString output  [10] = {"HIST"   , "HIST"   , "HIST"   , "HIST"   , "HIST"  , "HIST"  , "HIST"  ,  "HIST"  , "HIST"  ,  "HIST"  };
   Char_t  charge1 [10] = {'+'      , '-'      , '+'      , '-'      , '+'     , '-'     , '+'     ,  '-'     , '+'     ,  '-'     };
   Char_t  charge2 [10] = {'-'      , '+'      , '-'      , '+'      , '+'     , '-'     , '-'     ,  '+'     , '-'     ,  '+'     };
   Int_t   cutID1  [10] = { iCutK   ,  iCutK   ,  iCutK   ,  iCutK   ,  iCutK  ,  iCutK  ,  iCutK  ,   iCutK  ,  iCutK  ,   iCutK  };
   Int_t   cutID2  [10] = { iCutPi  ,  iCutPi  ,  iCutPi  ,  iCutPi  ,  iCutPi ,  iCutPi ,  iCutPi ,   iCutPi ,  iCutPi ,   iCutPi };
   
   for (Int_t i = 0; i < 10; i++) {
      if (!use[i]) continue;
      // create output
      AliRsnMiniOutput *out = task->CreateOutput(Form("kstar_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      out->SetCutID(0, cutID1[i]);
      out->SetCutID(1, cutID2[i]);
      out->SetDaughter(0, AliRsnDaughter::kKaon);
      out->SetDaughter(1, AliRsnDaughter::kPion);
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetMotherPDG(313);
      out->SetMotherMass(0.896);
      // pair cuts
      out->SetPairCuts(cutsPair);
      // axis X: invmass (or resolution)
      if (useIM[i]) 
         out->AddAxis(imID, 90, 0.6, 1.5);
      else
         out->AddAxis(resID, 200, -0.02, 0.02);
      // axis Y: transverse momentum
      out->AddAxis(ptID, 100, 0.0, 10.0);
   }
   
   return kTRUE;
}
