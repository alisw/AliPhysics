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
Bool_t ConfigSigmaStar
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
   
   // quality cuts
   AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("qualityDaughterLambda");
   
   esdTrackCuts->SetAcceptKinkDaughters(0); // 0 = kFalse
   //esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
   esdTrackCuts->SetMaxChi2PerClusterTPC(4);
   esdTrackCuts->SetMinNClustersTPC(70);
   esdTrackCuts->SetRequireTPCRefit();
   
   // cut lambda
   AliRsnCutV0 *cutLambda = new AliRsnCutV0("cutLambda", kLambda0);
   cutLambda->SetESDtrackCuts(esdTrackCuts);
   // cut set
   AliRsnCutSet *cutSetLambda = new AliRsnCutSet("setLambda", AliRsnTarget::kDaughter);
   cutSetLambda->AddCut(cutLambda);
   cutSetLambda->SetCutScheme(cutLambda->GetName());
   // add to task
   Int_t iCutLambda = task->AddTrackCuts(cutSetLambda);
   
   // cut anti-AntiLambda
   AliRsnCutV0 *cutAntiLambda = new AliRsnCutV0("cutAntiLambda", kLambda0Bar);
   cutAntiLambda->SetESDtrackCuts(esdTrackCuts);
   // cut set
   AliRsnCutSet *cutSetAntiLambda = new AliRsnCutSet("setAntiLambda", AliRsnTarget::kDaughter);
   cutSetAntiLambda->AddCut(cutAntiLambda);
   cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
   // add to task
   Int_t iCutAntiLambda = task->AddTrackCuts(cutSetAntiLambda);
   
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
   Bool_t   use     [12] = { 1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             };
   Bool_t   useIM   [12] = { 1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             };
   TString  name    [12] = {"SigmaP"   , "SigmaM"   , "ASigmaP"      , "ASigmaM"      , "SigmaPmix", "SigmaMmix", "ASigmaPmix"   , "ASigmaMmix"   , "SigmaPt"  , "SigmaMt"  , "ASigmaPt"     , "ASigmaMt"     };
   TString  comp    [12] = {"PAIR"     , "PAIR"     , "PAIR"         , "PAIR"         , "MIX"      , "MIX"      , "MIX"          , "MIX"          , "TRUE"     , "TRUE"     , "TRUE"         , "TRUE"         };
   TString  output  [12] = {"HIST"     , "HIST"     , "HIST"         , "HIST"         , "HIST"     , "HIST"     , "HIST"         , "HIST"         , "HIST"     , "HIST"     , "HIST"         , "HIST"         };
   Char_t   charge1 [12] = {'0'        , '0'        , '0'            , '0'            , '0'        , '0'        , '0'            , '0'            , '0'        , '0'        , '0'            , '0'            };
   Char_t   charge2 [12] = {'+'        , '-'        , '-'            , '+'            , '+'        , '-'        , '-'            , '+'            , '+'        , '-'        , '-'            , '+'            };
   Int_t    cutID1  [12] = { iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda,  iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda,  iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda};
   Int_t    cutID2  [12] = { iCutPi    ,  iCutPi    ,  iCutPi        ,  iCutPi        ,  iCutPi    ,  iCutPi    ,  iCutPi        ,  iCutPi        ,  iCutPi    ,  iCutPi    ,  iCutPi        ,  iCutPi        };
   Int_t    ipdg    [12] = { 3224      ,  3114      , -3224          , -3114          ,  3224      ,  3114      , -3224          , -3114          ,  3224      ,  3114      , -3224          , -3114          };
   Double_t mass    [12] = { 1382.3    ,  1387.4    ,  1382.3        ,  1387.4        ,  1382.3    ,  1387.4    ,  1382.3        ,  1387.4        ,  1382.3    ,  1387.4    ,  1382.3        ,  1387.4        };
   
   for (Int_t i = 0; i < 12; i++) {
      if (!use[i]) continue;
      // create output
      AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastar_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      out->SetCutID(0, cutID1[i]);
      out->SetCutID(1, cutID2[i]);
      out->SetDaughter(0, AliRsnDaughter::kLambda);
      out->SetDaughter(1, AliRsnDaughter::kPion);
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetMotherPDG(ipdg[i]);
      out->SetMotherMass(mass[i]);
      // pair cuts
      out->SetPairCuts(cutsPair);
      // axis X: invmass (or resolution)
      if (useIM[i]) 
         out->AddAxis(imID, 800, 1.2, 2.0);
      else
         out->AddAxis(resID, 200, -0.02, 0.02);
      // axis Y: transverse momentum
      out->AddAxis(ptID, 100, 0.0, 10.0);
   }
   
   return kTRUE;
}
