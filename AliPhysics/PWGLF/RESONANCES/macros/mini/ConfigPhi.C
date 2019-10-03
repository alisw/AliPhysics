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
Bool_t ConfigPhi
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
   
   // integrated kaon cut
   AliRsnCutKaonForPhi2010PP *cutStd = new AliRsnCutKaonForPhi2010PP("cutStdPP");
   // cut set
   AliRsnCutSet *cutSetStd = new AliRsnCutSet("kaonForPhi", AliRsnTarget::kDaughter);
   cutSetStd->AddCut(cutStd);
   cutSetStd->SetCutScheme(cutStd->GetName());
   // add to task
   Int_t icut = task->AddTrackCuts(cutSetStd);
   
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
   Bool_t  use     [6] = { 1      ,  1      ,  1      ,  1      ,  isMC   ,  isMC     };
   Bool_t  useIM   [6] = { 1      ,  1      ,  1      ,  1      ,  1      ,  0        };
   TString name    [6] = {"Unlike", "Mixing", "LikePP", "LikeMM", "Trues" , "Res"     };
   TString comp    [6] = {"PAIR"  , "MIX"   , "PAIR"  , "PAIR"  , "TRUE"  , "TRUE"    };
   TString output  [6] = {"HIST"  , "HIST"  , "HIST"  , "HIST"  , "HIST"  , "HIST"    };
   Char_t  charge1 [6] = {'+'     , '+'     , '+'     , '-'     , '+'     , '+'       };
   Char_t  charge2 [6] = {'-'     , '-'     , '+'     , '-'     , '-'     , '-'       };
   Int_t   cutID   [6] = { icut   ,  icut   ,  icut   ,  icut   ,  icut   ,  icut     };
   
   for (Int_t i = 0; i < 6; i++) {
      if (!use[i]) continue;
      // create output
      AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
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
      // axis X: invmass (or resolution)
      if (useIM[i]) 
         out->AddAxis(imID, 500, 0.9,  1.4);
      else
         out->AddAxis(resID, 200, -0.02, 0.02);
      // axis Y: transverse momentum
      out->AddAxis(ptID, 100, 0.0, 10.0);
   }
   
   return kTRUE;
}
