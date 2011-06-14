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
Bool_t ConfigPhi(AliRsnMiniAnalysisTask *task, Bool_t isMC, const char *suffix)
{
   // 
   // -- Define track cuts -------------------------------------------------------------------------
   //
   
   // integrated kaon cut
   AliRsnCutKaonForPhi2010PP *cutStd = new AliRsnCutKaonForPhi2010PP("cutStd");
   // cut set
   AliRsnCutSet *cutSetStd = new AliRsnCutSet("kaonForPhi", AliRsnTarget::kDaughter);
   cutSetStd->AddCut(cutStd);
   cutSetStd->SetCutScheme(cutStd->GetName());
   // add to task
   Int_t iCutStd = task->AddTrackCuts(cutSetStd);
   
   //
   // -- Values ------------------------------------------------------------------------------------
   //
   
   // invariant mass
   Int_t imID = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
   
   // transverse momentum
   Int_t ptID = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
   
   //
   // -- Pair cuts ---------------------------------------------------------------------------------
   //
   
   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   cutY->SetRangeD(-0.5, 0.5);
   
   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   
   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //
   
   // use an array for more compact writing, which are different on mixing and charges
   // [0] = unlike
   // [1] = mixing
   // [2] = like ++
   // [3] = like --
   Bool_t  use     [5] = { 1      ,  1      ,  1      ,  1      ,  1        };
   TString name    [5] = {"Unlike", "Mixing", "LikePP", "LikeMM", "Trues"   };
   TString comp    [5] = {"PAIR"  , "MIX"   , "PAIR"  , "PAIR"  , "TRUE"    };
   TString output  [5] = {"HIST"  , "HIST"  , "HIST"  , "HIST"  , "HIST"    };
   Char_t  charge1 [5] = {'+'     , '+'     , '+'     , '-'     , '+'       };
   Char_t  charge2 [5] = {'-'     , '-'     , '+'     , '-'     , '-'       };
   Int_t   cutID   [5] = { iCutStd,  iCutStd,  iCutStd,  iCutStd,  iCutStd  };
   
   for (Int_t i = 0; i < 5; i++) {
      if (!use[i]) continue;
      // create output
      AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s_%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      out->SetCutID(0, cutID[i]);
      out->SetCutID(1, cutID[i]);
      out->SetDaughter(0, AliRsnDaughter::kKaon);
      out->SetDaughter(1, AliRsnDaughter::kKaon);
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetMotherPDG(333);
      out->SetMotherMass(1.019455);
      // pair cuts
      out->SetPairCuts(cutsPair);
      // binnings
      out->AddAxis(imID, 500, 0.9,  1.4);
      out->AddAxis(ptID, 200, 0.0, 10.0);
   }
   
   return kTRUE;
}
