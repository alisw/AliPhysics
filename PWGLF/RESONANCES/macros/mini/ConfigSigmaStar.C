//
// *** Configuration script for Sigma*->Lambda-Pi analysis with 2010 runs ***
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
   Bool_t		   isPP, 
   Bool_t                  isMC,
   Float_t                 piPIDCut,
   Float_t                 pPIDCut,
   Int_t                  aodFilterBit,
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
   AliRsnCutDaughterSigmaStar2010PP *cutPi = new AliRsnCutDaughterSigmaStar2010PP("cutPionForSigmaStar", AliPID::kPion);
   cutPi->SetPIDCut(piPIDCut);
   AliRsnCutTrackQuality *cutQuality = (AliRsnCutTrackQuality*) cutPi->CutQuality();
   cutQuality->SetDCARmax(0.05);	         
    
   // cut set
   AliRsnCutSet *cutSetPi = new AliRsnCutSet("setPionForSigmaStar", AliRsnTarget::kDaughter);
   cutSetPi->AddCut(cutPi);
   cutSetPi->SetCutScheme(cutPi->GetName());
   // add to task
   Int_t iCutPi = task->AddTrackCuts(cutSetPi);
   
   // added by EF
   AliRsnCutTrackQuality *cq = (AliRsnCutTrackQuality*) cutPi->CutQuality();
   cq->Print();
   cq->SetAODTestFilterBit(-1);
   
   // quality cuts
   AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("qualityDaughterLambda");
   
   esdTrackCuts->SetAcceptKinkDaughters(0); // 0 = kFalse
   //esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
   esdTrackCuts->SetMaxChi2PerClusterTPC(4);
   esdTrackCuts->SetMinNClustersTPC(70);
   esdTrackCuts->SetRequireTPCRefit();
   
   // cut lambda
   //AliRsnCutV0 *cutLambda = new AliRsnCutV0("cutLambda", kLambda0);
   AliRsnCutV0 *cutLambda = new AliRsnCutV0("cutLambda", kLambda0, AliPID::kProton, AliPID::kPion);
   cutLambda->SetESDtrackCuts(esdTrackCuts);
   cutLambda->SetTolerance(0.01);
   cutLambda->SetMaxDCAVertex(0.3);
   cutLambda->SetMinCosPointingAngle(0.99);
   cutLambda->SetMaxDaughtersDCA(0.5);
   cutLambda->SetMaxRapidity(0.5);
   cutLambda->SetAODTestFilterBit(aodFilterBit);
   cutLambda->SetPIDCut1(pPIDCut);
   cutLambda->SetPIDCut2(piPIDCut);
   cutLambda->SetPIDCut3(piPIDCut);
   
   // cut set
   AliRsnCutSet *cutSetLambda = new AliRsnCutSet("setLambda", AliRsnTarget::kDaughter);
   cutSetLambda->AddCut(cutLambda);
   cutSetLambda->SetCutScheme(cutLambda->GetName());
   
   // add to task
   Int_t iCutLambda = task->AddTrackCuts(cutSetLambda);
   
   // cut anti-AntiLambda
   //AliRsnCutV0 *cutAntiLambda = new AliRsnCutV0("cutAntiLambda", kLambda0Bar);
   AliRsnCutV0 *cutAntiLambda = new AliRsnCutV0("cutAntiLambda", kLambda0Bar, AliPID::kProton, AliPID::kPion);
   cutAntiLambda->SetESDtrackCuts(esdTrackCuts);
   cutAntiLambda->SetTolerance(0.01);
   cutAntiLambda->SetMaxDCAVertex(0.3);
   cutAntiLambda->SetMinCosPointingAngle(0.99);
   cutAntiLambda->SetMaxDaughtersDCA(0.5);
   cutAntiLambda->SetMaxRapidity(0.5);
   cutAntiLambda->SetAODTestFilterBit(aodFilterBit);
   cutAntiLambda->SetPIDCut1(pPIDCut);
   cutAntiLambda->SetPIDCut2(piPIDCut);
   cutAntiLambda->SetPIDCut3(piPIDCut);
   
   // cut set
   AliRsnCutSet *cutSetAntiLambda = new AliRsnCutSet("setAntiLambda", AliRsnTarget::kDaughter);
   cutSetAntiLambda->AddCut(cutAntiLambda);
   cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
   // add to task
   Int_t iCutAntiLambda = task->AddTrackCuts(cutSetAntiLambda);
   
   //######################################################################################################
   //Loose Cuts
   
   // integrated pion cut
   AliRsnCutDaughterSigmaStar2010PP *cutPiLoose = new AliRsnCutDaughterSigmaStar2010PP("cutPionForSigmaStarLoose", AliPID::kPion);
   cutPiLoose->SetPIDCut(piPIDCut);
   AliRsnCutTrackQuality *cutQuality = (AliRsnCutTrackQuality*) cutPi->CutQuality();
   cutQuality->SetDCARmax(0.072);
   
   // cut set
   AliRsnCutSet *cutSetPiLoose = new AliRsnCutSet("setPionForSigmaStarLoose", AliRsnTarget::kDaughter);
   cutSetPiLoose->AddCut(cutPiLoose);
   cutSetPiLoose->SetCutScheme(cutPiLoose->GetName());
   // add to task
   Int_t iCutPiLoose = task->AddTrackCuts(cutSetPiLoose);
   
   
   // cut lambda
   //AliRsnCutV0 *cutLambdaLoose = new AliRsnCutV0("cutLambdaLoose", kLambda0);
   AliRsnCutV0 *cutLambdaLoose = new AliRsnCutV0("cutLambdaLoose", kLambda0, AliPID::kProton, AliPID::kPion);
   cutLambdaLoose->SetESDtrackCuts(esdTrackCuts);
   cutLambdaLoose->SetTolerance(0.013);
   cutLambdaLoose->SetMaxDCAVertex(0.4);
   cutLambdaLoose->SetMinCosPointingAngle(0.82);
   cutLambdaLoose->SetMaxDaughtersDCA(0.67);
   cutLambdaLoose->SetMaxRapidity(0.5);
   cutLambdaLoose->SetAODTestFilterBit(aodFilterBit);
   cutLambdaLoose->SetPIDCut1(pPIDCut);
   cutLambdaLoose->SetPIDCut2(piPIDCut);
   cutLambdaLoose->SetPIDCut3(piPIDCut);
   
   // cut set
   AliRsnCutSet *cutSetLambdaLoose = new AliRsnCutSet("setLambdaLoose", AliRsnTarget::kDaughter);
   cutSetLambdaLoose->AddCut(cutLambdaLoose);
   cutSetLambdaLoose->SetCutScheme(cutLambdaLoose->GetName());
   
   // add to task
   Int_t iCutLambdaLoose = task->AddTrackCuts(cutSetLambdaLoose);
   
   // cut anti-AntiLambda
   //AliRsnCutV0 *cutAntiLambdaLoose = new AliRsnCutV0("cutAntiLambdaLoose", kLambda0Bar);
   AliRsnCutV0 *cutAntiLambdaLoose = new AliRsnCutV0("cutAntiLambdaLoose", kLambda0Bar, AliPID::kProton, AliPID::kPion);
   cutAntiLambdaLoose->SetESDtrackCuts(esdTrackCuts);
   cutAntiLambdaLoose->SetTolerance(0.013);
   cutAntiLambdaLoose->SetMaxDCAVertex(0.4);
   cutAntiLambdaLoose->SetMinCosPointingAngle(0.982);
   cutAntiLambdaLoose->SetMaxDaughtersDCA(0.67);
   cutAntiLambdaLoose->SetMaxRapidity(0.5);
   cutAntiLambdaLoose->SetAODTestFilterBit(aodFilterBit);
   cutAntiLambdaLoose->SetPIDCut1(pPIDCut);
   cutAntiLambdaLoose->SetPIDCut2(piPIDCut);
   cutAntiLambdaLoose->SetPIDCut3(piPIDCut);
   
   // cut set
   AliRsnCutSet *cutSetAntiLambdaLoose = new AliRsnCutSet("setAntiLambdaLoose", AliRsnTarget::kDaughter);
   cutSetAntiLambdaLoose->AddCut(cutAntiLambdaLoose);
   cutSetAntiLambdaLoose->SetCutScheme(cutAntiLambdaLoose->GetName());
   // add to task
   Int_t iCutAntiLambdaLoose = task->AddTrackCuts(cutSetAntiLambdaLoose);
   

   //######################################################################################################  
   
   
   //######################################################################################################
   //Tight Cuts
   
   // integrated pion cut
   AliRsnCutDaughterSigmaStar2010PP *cutPiTight = new AliRsnCutDaughterSigmaStar2010PP("cutPionForSigmaStarTight", AliPID::kPion);
   cutPiTight->SetPIDCut(piPIDCut);
   AliRsnCutTrackQuality *cutQuality = (AliRsnCutTrackQuality*) cutPi->CutQuality();
   cutQuality->SetDCARmax(0.036);
   
   // cut set
   AliRsnCutSet *cutSetPiTight = new AliRsnCutSet("setPionForSigmaStarTight", AliRsnTarget::kDaughter);
   cutSetPiTight->AddCut(cutPiTight);
   cutSetPiTight->SetCutScheme(cutPiTight->GetName());
   // add to task
   Int_t iCutPiTight = task->AddTrackCuts(cutSetPiTight);
   
   
   // cut lambda
   //AliRsnCutV0 *cutLambdaTight = new AliRsnCutV0("cutLambdaTight", kLambda0);
   AliRsnCutV0 *cutLambdaTight = new AliRsnCutV0("cutLambdaTight", kLambda0, AliPID::kProton, AliPID::kPion);
   cutLambdaTight->SetESDtrackCuts(esdTrackCuts);
   cutLambdaTight->SetTolerance(0.0076);
   cutLambdaTight->SetMaxDCAVertex(0.22);
   cutLambdaTight->SetMinCosPointingAngle(0.995);
   cutLambdaTight->SetMaxDaughtersDCA(0.35);
   cutLambdaTight->SetMaxRapidity(0.5);
   cutLambdaTight->SetAODTestFilterBit(aodFilterBit);
   cutLambdaTight->SetPIDCut1(pPIDCut);
   cutLambdaTight->SetPIDCut2(piPIDCut);
   cutLambdaTight->SetPIDCut3(piPIDCut);
   
   // cut set
   AliRsnCutSet *cutSetLambdaTight = new AliRsnCutSet("setLambdaTight", AliRsnTarget::kDaughter);
   cutSetLambdaTight->AddCut(cutLambdaTight);
   cutSetLambdaTight->SetCutScheme(cutLambdaTight->GetName());
   
   // add to task
   Int_t iCutLambdaTight = task->AddTrackCuts(cutSetLambdaTight);
   
   // cut anti-AntiLambda
   //AliRsnCutV0 *cutAntiLambdaTight = new AliRsnCutV0("cutAntiLambdaTight", kLambda0Bar);
   AliRsnCutV0 *cutAntiLambdaTight = new AliRsnCutV0("cutAntiLambdaTight", kLambda0Bar, AliPID::kProton, AliPID::kPion);
   cutAntiLambdaTight->SetESDtrackCuts(esdTrackCuts);
   cutAntiLambdaTight->SetTolerance(0.076);
   cutAntiLambdaTight->SetMaxDCAVertex(0.22);
   cutAntiLambdaTight->SetMinCosPointingAngle(0.995);
   cutAntiLambdaTight->SetMaxDaughtersDCA(0.35);
   cutAntiLambdaTight->SetMaxRapidity(0.5);
   cutAntiLambdaTight->SetAODTestFilterBit(aodFilterBit);
   cutAntiLambdaTight->SetPIDCut1(pPIDCut);
   cutAntiLambdaTight->SetPIDCut2(piPIDCut);
   cutAntiLambdaTight->SetPIDCut3(piPIDCut);
   
   // cut set
   AliRsnCutSet *cutSetAntiLambdaTight = new AliRsnCutSet("setAntiLambdaTight", AliRsnTarget::kDaughter);
   cutSetAntiLambdaTight->AddCut(cutAntiLambdaTight);
   cutSetAntiLambdaTight->SetCutScheme(cutAntiLambdaTight->GetName());
   // add to task
   Int_t iCutAntiLambdaTight = task->AddTrackCuts(cutSetAntiLambdaTight);
   

   
   //######################################################################################################  
    
   //
   // -- Values ------------------------------------------------------------------------------------
   //
   
   /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
   /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
   /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
   /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
   /* rapidity         */ Int_t rapID  = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
   
   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //
   
   // use an array for more compact writing, which are different on mixing and charges
   // [0] = unlike
   // [1] = mixing
   // [2] = like ++
   // [3] = like --
   Bool_t   use     [18] = { 1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1          ,  1              ,  1              ,  1              ,  1              ,  1              };
   Bool_t   useIM   [18] = { 1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1          ,  1              ,  1              ,  1              ,  1              ,  1              };
   TString  name    [18] = {"SigmaP"   , "SigmaM"   , "ASigmaP"      , "ASigmaM"      , "SigmaPmix", "SigmaMmix", "ASigmaPmix"   , "ASigmaMmix"   , "SigmaPt"  , "SigmaMt"  , "ASigmaPt"     , "ASigmaMt"     , "XiM"       , "XiP"           , "Lambda1520P"   , "Lambda1520M"   , "Lambda1520PBar", "Lambda1520MBar"};
   TString  comp    [18] = {"PAIR"     , "PAIR"     , "PAIR"         , "PAIR"         , "MIX"      , "MIX"      , "MIX"          , "MIX"          , "TRUE"     , "TRUE"     , "TRUE"         , "TRUE"         , "TRUE"      , "TRUE"          , "TRUE"          , "TRUE"          , "TRUE"          , "TRUE"          };
   TString  output  [18] = {"HIST"     , "HIST"     , "HIST"         , "HIST"         , "HIST"     , "HIST"     , "HIST"         , "HIST"         , "HIST"     , "HIST"     , "HIST"         , "HIST"         , "HIST"      , "HIST"          , "HIST"          , "HIST"          , "HIST"          , "HIST"          };
   Char_t   charge1 [18] = {'0'        , '0'        , '0'            , '0'            , '0'        , '0'        , '0'            , '0'            , '0'        , '0'        , '0'            , '0'            , '0'         , '0'             , '0'             , '0'             , '0'             , '0'             };
   Char_t   charge2 [18] = {'+'        , '-'        , '-'            , '+'            , '+'        , '-'        , '-'            , '+'            , '+'        , '-'        , '-'            , '+'            , '-'         , '+'             , '+'             , '-'             , '-'             , '+'             };
   Int_t    cutID1  [18] = { iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda,  iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda,  iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda,  iCutLambda ,  iCutAntiLambda ,  iCutLambda     ,  iCutLambda     ,  iCutAntiLambda ,  iCutAntiLambda };
   Int_t    cutID2  [18] = { iCutPi    ,  iCutPi    ,  iCutPi        ,  iCutPi        ,  iCutPi    ,  iCutPi    ,  iCutPi        ,  iCutPi        ,  iCutPi    ,  iCutPi    ,  iCutPi        ,  iCutPi        ,  iCutPi     ,  iCutPi         ,  iCutPi         ,  iCutPi         ,  iCutPi         ,  iCutPi         };
   Int_t    ipdg    [18] = { 3224      ,  3114      , -3224          , -3114          ,  3224      ,  3114      , -3224          , -3114          ,  3224      ,  3114      , -3224          , -3114          ,  3312       , -3312           ,  3124           ,  3124           , -3124           , -3124           };
   Double_t mass    [18] = { 1.3828    ,  1.3872    ,  1.3828        ,  1.3872        ,  1.3828    ,  1.3872    ,  1.3828        ,  1.3872        ,  1.3828    ,  1.3872    ,  1.3828        ,  1.3872        ,  1.32171    ,  1.32171        ,  1.5195         ,  1.5195         ,  1.5195         ,  1.5195         };
   
   for (Int_t i = 0; i < 18; i++) {
      if (!use[i]) continue;
      if (!isPP) output[i] = "SPARSE";
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
      // axis X: invmass
      if (useIM[i]) 
         out->AddAxis(imID, 800, 1.2, 2.0);
      // axis Y: transverse momentum
         out->AddAxis(ptID, 100, 0.0, 10.0);
      // axis Z: rapidity
         //out->AddAxis(rapID, 160, -0.8, 0.8);
	 
      if (!isPP) out->AddAxis(centID, 100, 0.0, 100.0);
      
    } 
   
   
   //############################################################################
   //Loose Cuts
   
   TString  nameLoose    [18] = {"SigmaPLoose"   , "SigmaMLoose"   , "ASigmaPLoose"      , "ASigmaMLoose"      , "SigmaPmixLoose", "SigmaMmixLoose", "ASigmaPmixLoose"   , "ASigmaMmixLoose"   , "SigmaPtLoose"  , "SigmaMtLoose"  , "ASigmaPtLoose"     , "ASigmaMtLoose"     , "XiMLoose"        , "XiPLoose"            , "Lambda1520PLoose" , "Lambda1520MLoose", "Lambda1520PBarLoose", "Lambda1520MBarLoose"};
   Int_t    cutID1Loose  [18] = { iCutLambdaLoose,  iCutLambdaLoose,  iCutAntiLambdaLoose,  iCutAntiLambdaLoose,  iCutLambdaLoose,  iCutLambdaLoose,  iCutAntiLambdaLoose,  iCutAntiLambdaLoose,  iCutLambdaLoose,  iCutLambdaLoose,  iCutAntiLambdaLoose,  iCutAntiLambdaLoose,  iCutLambdaLoose  ,  iCutAntiLambdaLoose  ,  iCutLambdaLoose   ,  iCutLambdaLoose  ,  iCutAntiLambdaLoose ,  iCutAntiLambdaLoose };
   Int_t    cutID2Loose  [18] = { iCutPiLoose    ,  iCutPiLoose    ,  iCutPiLoose        ,  iCutPiLoose        ,  iCutPiLoose    ,  iCutPiLoose    ,  iCutPiLoose        ,  iCutPiLoose        ,  iCutPiLoose    ,  iCutPiLoose    ,  iCutPiLoose        ,  iCutPiLoose        ,  iCutPiLoose      ,  iCutPiLoose          ,  iCutPiLoose       ,  iCutPiLoose      ,  iCutPiLoose         ,  iCutPiLoose         };
   
   
   for (Int_t i = 0; i < 18; i++) {
      if (!use[i]) continue;
      if (!isPP) output[i] = "SPARSE";
      // create output
      AliRsnMiniOutput *outLoose = task->CreateOutput(Form("sigmastar_%s%s", nameLoose[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      outLoose->SetCutID(0, cutID1Loose[i]);
      outLoose->SetCutID(1, cutID2Loose[i]);
      outLoose->SetDaughter(0, AliRsnDaughter::kLambda);
      outLoose->SetDaughter(1, AliRsnDaughter::kPion);
      outLoose->SetCharge(0, charge1[i]);
      outLoose->SetCharge(1, charge2[i]);
      outLoose->SetMotherPDG(ipdg[i]);
      outLoose->SetMotherMass(mass[i]);
      // pair cuts
      out->SetPairCuts(cutsPair);
      // axis X: invmass
      if (useIM[i]) 
         outLoose->AddAxis(imID, 800, 1.2, 2.0);
      // axis Y: transverse momentum
         outLoose->AddAxis(ptID, 100, 0.0, 10.0);
      // axis Z: rapidity
         //out->AddAxis(rapID, 160, -0.8, 0.8);
      if (!isPP) outLoose->AddAxis(centID, 100, 0.0, 100.0);
       
   }
   
   
   //############################################################################
   
   
   //############################################################################
   //Tight Cuts
   
   TString  nameTight    [18] = {"SigmaPTight"   , "SigmaMTight"   , "ASigmaPTight"      , "ASigmaMTight"      , "SigmaPmixTight", "SigmaMmixTight", "ASigmaPmixTight"   , "ASigmaMmixTight"   , "SigmaPtTight"  , "SigmaMtTight"  , "ASigmaPtTight"     , "ASigmaMtTight"     , "XiMTight"        , "XiPTight"            , "Lambda1520PTight" , "Lambda1520MTight", "Lambda1520PBarTight", "Lambda1520MBarTight"};
   Int_t    cutID1Tight  [18] = { iCutLambdaTight,  iCutLambdaTight,  iCutAntiLambdaTight,  iCutAntiLambdaTight,  iCutLambdaTight,  iCutLambdaTight,  iCutAntiLambdaTight,  iCutAntiLambdaTight,  iCutLambdaTight,  iCutLambdaTight,  iCutAntiLambdaTight,  iCutAntiLambdaTight,  iCutLambdaTight  ,  iCutAntiLambdaTight  ,  iCutLambdaTight   ,  iCutLambdaTight  ,  iCutAntiLambdaTight ,  iCutAntiLambdaTight };
   Int_t    cutID2Tight  [18] = { iCutPiTight    ,  iCutPiTight    ,  iCutPiTight        ,  iCutPiTight        ,  iCutPiTight    ,  iCutPiTight    ,  iCutPiTight        ,  iCutPiTight        ,  iCutPiTight    ,  iCutPiTight    ,  iCutPiTight        ,  iCutPiTight        ,  iCutPiTight      ,  iCutPiTight          ,  iCutPiTight       ,  iCutPiTight      ,  iCutPiTight         ,  iCutPiTight         };
   
   
   for (Int_t i = 0; i < 18; i++) {
      if (!use[i]) continue;
      if (!isPP) output[i] = "SPARSE";
      // create output
      AliRsnMiniOutput *outTight = task->CreateOutput(Form("sigmastar_%s%s", nameTight[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      outTight->SetCutID(0, cutID1Tight[i]);
      outTight->SetCutID(1, cutID2Tight[i]);
      outTight->SetDaughter(0, AliRsnDaughter::kLambda);
      outTight->SetDaughter(1, AliRsnDaughter::kPion);
      outTight->SetCharge(0, charge1[i]);
      outTight->SetCharge(1, charge2[i]);
      outTight->SetMotherPDG(ipdg[i]);
      outTight->SetMotherMass(mass[i]);
      // pair cuts
      outTight->SetPairCuts(cutsPair);
      // axis X: invmass
      if (useIM[i]) 
         outTight->AddAxis(imID, 800, 1.2, 2.0);
      // axis Y: transverse momentum
         outTight->AddAxis(ptID, 100, 0.0, 10.0);
      // axis Z: rapidity
         //out->AddAxis(rapID, 160, -0.8, 0.8);
      if (!isPP) outTight->AddAxis(centID, 100, 0.0, 100.0);
       
   }
   
   
   AddMonitorOutput(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput(cutSetAntiLambda->GetMonitorOutput());
   
   
   //############################################################################
   
   if (isMC) {
   
   TString mode = "HIST";
   if (!isPP) mode = "SPARSE";
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarP_TrueMC%s", suffix), mode.Data(), "MOTHER");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kLambda);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetMotherPDG(3224);
   out->SetMotherMass(1.3828);
   // pair cuts
   out->SetPairCuts(cutsPair);
   // binnings
   out->AddAxis(imID, 800, 1.2, 2.0);
   out->AddAxis(ptID, 100, 0.0, 10.0);
   //out->AddAxis(rapID, 160, -0.8, 0.8);
   if (!isPP) out->AddAxis(centID, 100, 0.0, 100.0);
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarM_TrueMC%s", suffix), mode.Data(), "MOTHER");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kLambda);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetMotherPDG(3114);
   out->SetMotherMass(1.3872);
   // pair cuts
   out->SetPairCuts(cutsPair);
   // binnings
   out->AddAxis(imID, 800, 1.2, 2.0);
   out->AddAxis(ptID, 100, 0.0, 10.0);
   //out->AddAxis(rapID, 160, -0.8, 0.8);
   if (!isPP) out->AddAxis(centID, 100, 0.0, 100.0);
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarPBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kLambda);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetMotherPDG(-3224);
   out->SetMotherMass(1.3828);
   // pair cuts
   out->SetPairCuts(cutsPair);
   // binnings
   out->AddAxis(imID, 800, 1.2, 2.0);
   out->AddAxis(ptID, 100, 0.0, 10.0);
   //out->AddAxis(rapID, 160, -0.8, 0.8);
   if (!isPP) out->AddAxis(centID, 100, 0.0, 100.0);
   
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarMBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kLambda);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetMotherPDG(-3114);
   out->SetMotherMass(1.3872);
   // pair cuts
   out->SetPairCuts(cutsPair);
   // binnings
   out->AddAxis(imID, 800, 1.2, 2.0);
   out->AddAxis(ptID, 100, 0.0, 10.0);
   //out->AddAxis(rapID, 160, -0.8, 0.8);
   if (!isPP) out->AddAxis(centID, 100, 0.0, 100.0);
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput(Form("XiP_TrueMC%s", suffix), mode.Data(), "MOTHER");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kLambda);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetMotherPDG(-3312);
   out->SetMotherMass(1.32171);
   // pair cuts
   out->SetPairCuts(cutsPair);
   // binnings
   out->AddAxis(imID, 800, 1.2, 2.0);
   out->AddAxis(ptID, 100, 0.0, 10.0);
   //out->AddAxis(rapID, 160, -0.8, 0.8);
   if (!isPP) out->AddAxis(centID, 100, 0.0, 100.0);
   
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput(Form("XiM_TrueMC%s", suffix), mode.Data(), "MOTHER");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kLambda);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetMotherPDG(3312);
   out->SetMotherMass(1.32171);
   // pair cuts
   out->SetPairCuts(cutsPair);
   // binnings
   out->AddAxis(imID, 800, 1.2, 2.0);
   out->AddAxis(ptID, 100, 0.0, 10.0);
   //out->AddAxis(rapID, 160, -0.8, 0.8);
   if (!isPP) out->AddAxis(centID, 100, 0.0, 100.0);
   
   
   AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520P_TrueMC%s", suffix), mode.Data(), "MOTHER");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kLambda);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetCharge(0, 0);
   out->SetCharge(1, 1);
   out->SetMotherPDG(3124);
   out->SetMotherMass(1.5195);
   // pair cuts
   out->SetPairCuts(cutsPair);
   // binnings
   out->AddAxis(imID, 800, 1.2, 2.0);
   out->AddAxis(ptID, 100, 0.0, 10.0);
   //out->AddAxis(rapID, 160, -0.8, 0.8);
   if (!isPP) out->AddAxis(centID, 100, 0.0, 100.0);
   
   AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520M_TrueMC%s", suffix), mode.Data(), "MOTHER");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kLambda);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetCharge(0, 0);
   out->SetCharge(1, -1);
   out->SetMotherPDG(3124);
   out->SetMotherMass(1.5195);
   // pair cuts
   out->SetPairCuts(cutsPair);
   // binnings
   out->AddAxis(imID, 800, 1.2, 2.0);
   out->AddAxis(ptID, 100, 0.0, 10.0);
   //out->AddAxis(rapID, 160, -0.8, 0.8);
   if (!isPP) out->AddAxis(centID, 100, 0.0, 100.0);
   
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520PBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kLambda);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetCharge(0, 0);
   out->SetCharge(1, 1);
   out->SetMotherPDG(-3124);
   out->SetMotherMass(1.5195);
   // pair cuts
   out->SetPairCuts(cutsPair);
   // binnings
   out->AddAxis(imID, 800, 1.2, 2.0);
   out->AddAxis(ptID, 100, 0.0, 10.0);
   //out->AddAxis(rapID, 160, -0.8, 0.8);
   if (!isPP) out->AddAxis(centID, 100, 0.0, 100.0);
   
   AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520MBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kLambda);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetCharge(0, 0);
   out->SetCharge(1, -1);
   out->SetMotherPDG(-3124);
   out->SetMotherMass(1.5195);
   // pair cuts
   out->SetPairCuts(cutsPair);
   // binnings
   out->AddAxis(imID, 800, 1.2, 2.0);
   out->AddAxis(ptID, 100, 0.0, 10.0);
   //out->AddAxis(rapID, 160, -0.8, 0.8);
   if (!isPP) out->AddAxis(centID, 100, 0.0, 100.0);
   
   
   
   }

   return kTRUE;
}



void AddMonitorOutput(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lm=0)
{

   // Mass
   AliRsnValueDaughter *axisMass = new AliRsnValueDaughter("m", AliRsnValueDaughter::kMass);
   axisMass->SetBins(0.7,1.5,0.05);

   // output: 2D histogram
   AliRsnListOutput *outMonitorM = new AliRsnListOutput("M", AliRsnListOutput::kHistoDefault);
   outMonitorM->AddValue(axisMass);

   // add outputs to loop
   if (mon) mon->Add(outMonitorM);
   if (lm) lm->AddOutput(outMonitorM);
  
}
