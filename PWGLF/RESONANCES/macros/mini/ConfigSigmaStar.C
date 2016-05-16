//
// *** Configuration script for Sigma*->Lambda-Pi analysis ***
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
   Int_t		   collSyst, 
   Bool_t                  isMC,
   Float_t                 piPIDCut,
   Float_t                 pPIDCut,
   Int_t                   aodFilterBit,
   Float_t                 trackDCAcut,
   Float_t                 massTol,
   Float_t                 lambdaDCA,
   Float_t                 lambdaCosPoinAn,
   Float_t                 lambdaDaughDCA,
   Int_t                   NTPCcluster,
   const char             *suffix,
   AliRsnCutSet           *cutsPair
)
{
   // manage suffix
   if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
 
   /////////////////////////////////////////////////////
   // selections for the pion from the decay of Sigma*
   /////////////////////////////////////////////////////
   //
   AliRsnCutDaughterSigmaStar2010PP *cutPi = new AliRsnCutDaughterSigmaStar2010PP("cutPionForSigmaStar", AliPID::kPion);
   cutPi->SetPIDCut(piPIDCut);    // fPIDCut used in IsSelected() after the call to cutQuality
   AliRsnCutTrackQuality *cutQuality = (AliRsnCutTrackQuality*) cutPi->CutQuality();
   //cutQuality->SetDefaults2011();
   cutQuality->SetDefaults2010(0,1);  // 1st par. not default (0 -> use TPC clusters). 2nd par. default (-> standard Pt and eta range)
   // SetDefaults2010 contains the following selections:
   //     SetPtRange(0.15, 1E+20);
   //     SetEtaRange(-0.8, 0.8);
   //     and from aliroot/master/src/ANALYSIS/ANALYSISalice/AliESDtrackCuts.cxx
   //     AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(1,0)
   //         esdTrackCuts->SetMinNClustersTPC(70);
   //         esdTrackCuts->SetMaxChi2PerClusterTPC(4);
   //         esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
   //         esdTrackCuts->SetRequireTPCRefit(kTRUE);
   //         esdTrackCuts->SetRequireITSRefit(kTRUE);
   //         esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
   //         esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");    // NB. With pt_min=0.15 (see above) -> DCAxy_max = 0.2560
   //         esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
   //         esdTrackCuts->SetMaxDCAToVertexZ(2);
   //         esdTrackCuts->SetDCAToVertex2D(kFALSE);
   //         esdTrackCuts->SetRequireSigmaToVertex(kFALSE);  
   //         esdTrackCuts->SetMaxChi2PerClusterITS(36);
   //  
   AliRsnCutSet *cutSetPi = new AliRsnCutSet("setPionForSigmaStar", AliRsnTarget::kDaughter);
   cutSetPi->AddCut(cutPi);
   cutSetPi->SetCutScheme(cutPi->GetName());
   Int_t iCutPi = task->AddTrackCuts(cutSetPi);
   //
   /////////////////////////////////////////////////////////////
   // selections for Lambda and for the daughters of Lambda 
   /////////////////////////////////////////////////////////////
   // 
   // selections for the proton and pion daugthers of Lambda and AntiLambda
   AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("qualityDaughterLambda");   
   esdTrackCuts->SetPtRange(0.15,1.E10);
   esdTrackCuts->SetEtaRange(-0.8,0.8);
   esdTrackCuts->SetRequireTPCRefit();
   esdTrackCuts->SetAcceptKinkDaughters(0); //
   esdTrackCuts->SetMinNClustersTPC(NTPCcluster);
   esdTrackCuts->SetMaxChi2PerClusterTPC(4);
   esdTrackCuts->SetMinDCAToVertexXY(0.15);   
   //
   /////////////////////////////////////////////////
   // selections for Lambda
   AliRsnCutV0 *cutLambda = new AliRsnCutV0("cutLambda", kLambda0, AliPID::kProton, AliPID::kPion);
   cutLambda->SetPIDCutProton(pPIDCut);       // PID for the proton daughter of Lambda
   cutLambda->SetPIDCutPion(piPIDCut);        // PID for the pion daughter of Lambda 
   cutLambda->SetESDtrackCuts(esdTrackCuts);  // all the other selections (defined above) for proton and pion daughters of Lambda
   cutLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
   cutLambda->SetMaxDCAVertex(lambdaDCA);
   cutLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
   cutLambda->SetTolerance(massTol);
   cutLambda->SetMaxRapidity(0.5);
   //
   AliRsnCutSet *cutSetLambda = new AliRsnCutSet("setLambda", AliRsnTarget::kDaughter);
   cutSetLambda->AddCut(cutLambda);
   cutSetLambda->SetCutScheme(cutLambda->GetName());
   Int_t iCutLambda = task->AddTrackCuts(cutSetLambda);
   //
   /////////////////////////////////////////////////
   // selections for AntiLambda
   AliRsnCutV0 *cutAntiLambda = new AliRsnCutV0("cutAntiLambda", kLambda0Bar, AliPID::kProton, AliPID::kPion);
   cutAntiLambda->SetPIDCutProton(pPIDCut);
   cutAntiLambda->SetPIDCutPion(piPIDCut);
   cutAntiLambda->SetESDtrackCuts(esdTrackCuts);
   cutAntiLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
   cutAntiLambda->SetMaxDCAVertex(lambdaDCA);
   cutAntiLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
   cutAntiLambda->SetTolerance(massTol);
   cutAntiLambda->SetMaxRapidity(0.5);
   // 
   AliRsnCutSet *cutSetAntiLambda = new AliRsnCutSet("setAntiLambda", AliRsnTarget::kDaughter);
   cutSetAntiLambda->AddCut(cutAntiLambda);
   cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
   Int_t iCutAntiLambda = task->AddTrackCuts(cutSetAntiLambda); 
   //
   /////////////////////////////////////////////////
   
   
   //######################################################################################################  
    
   //
   // -- Values ------------------------------------------------------------------------------------
   //
   
   /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
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
      if (collSyst) output[i] = "SPARSE";
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
      //  out->AddAxis(imID, 700, 1.2, 4.0);
      // axis Y: transverse momentum
	  out->AddAxis(ptID, 100, 0.0, 10.0);
	 //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
	 
      if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
      
    } 
    
   AddMonitorOutput_PionPt(cutSetPi->GetMonitorOutput());
   AddMonitorOutput_PionEta(cutSetPi->GetMonitorOutput());
   AddMonitorOutput_PionDCAxy(cutSetPi->GetMonitorOutput());
   AddMonitorOutput_PionDCAz(cutSetPi->GetMonitorOutput());
   AddMonitorOutput_PionPIDCut(cutSetPi->GetMonitorOutput());
   AddMonitorOutput_PionNTPC(cutSetPi->GetMonitorOutput());
   AddMonitorOutput_PionTPCchi2(cutSetPi->GetMonitorOutput());
   
   // AddMonitorOutput_LambdaP(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaPt(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaNegDaughPt(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaPosDaughPt(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaMass(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaDCA(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaRadius(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaDaughterDCA(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaCosPointAngle(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaProtonPID(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaPionPID(cutSetLambda->GetMonitorOutput());
   
   /*
   AddMonitorOutput_LambdaMass(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaP(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaPt(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaNegDaughPt(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaPosDaughPt(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaDCA(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaRadius(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaDaughterDCA(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaCosPointAngle(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaAntiProtonPID(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaAntiPionPID(cutSetAntiLambda->GetMonitorOutput());
   */

   if (isMC) {
     
     TString mode = "HIST";
     if (collSyst) mode = "SPARSE";
     
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
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
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
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
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
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
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
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
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
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
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
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
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
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
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
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
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
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
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
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     
   }
   
   return kTRUE;
}

void AddMonitorOutput_PionPt(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *ppt=0)
{

   // PionDCA
   AliRsnValueDaughter *axisPionPt = new AliRsnValueDaughter("pion_pt", AliRsnValueDaughter::kPt);
   axisPionPt->SetBins(0.,10.0,0.001);

   // output: 2D histogram
   AliRsnListOutput *outMonitorPionPt = new AliRsnListOutput("Pion_Pt", AliRsnListOutput::kHistoDefault);
   outMonitorPionPt->AddValue(axisPionPt);

   // add outputs to loop
   if (mon) mon->Add(outMonitorPionPt);
   if (ppt) ppt->AddOutput(outMonitorPionPt);
  
}

void AddMonitorOutput_PionEta(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *peta=0)
{

   // PionDCA
   AliRsnValueDaughter *axisPionEta = new AliRsnValueDaughter("pion_eta", AliRsnValueDaughter::kEta);
   axisPionEta->SetBins(-2.,2.,0.001);

   // output: 2D histogram
   AliRsnListOutput *outMonitorPionEta = new AliRsnListOutput("Pion_Eta", AliRsnListOutput::kHistoDefault);
   outMonitorPionEta->AddValue(axisPionEta);

   // add outputs to loop
   if (mon) mon->Add(outMonitorPionEta);
   if (peta) peta->AddOutput(outMonitorPionEta);
  
}

void AddMonitorOutput_PionDCAxy(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *pdcaxy=0)
{

   // PionDCA
   AliRsnValueDaughter *axisPionDCAxy = new AliRsnValueDaughter("pion_dcaxy", AliRsnValueDaughter::kDCAXY);
   axisPionDCAxy->SetBins(-0.5,0.5,0.001);

   // output: 2D histogram
   AliRsnListOutput *outMonitorPionDCAxy = new AliRsnListOutput("Pion_DCAxy", AliRsnListOutput::kHistoDefault);
   outMonitorPionDCAxy->AddValue(axisPionDCAxy);

   // add outputs to loop
   if (mon) mon->Add(outMonitorPionDCAxy);
   if (pdcaxy) pdcaxy->AddOutput(outMonitorPionDCAxy);
  
}

void AddMonitorOutput_PionDCAz(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *pdcaz=0)
{

   // PionDCA
   AliRsnValueDaughter *axisPionDCAz = new AliRsnValueDaughter("pion_dcaz", AliRsnValueDaughter::kDCAZ);
   axisPionDCAz->SetBins(-2.5,2.5,0.005);

   // output: 2D histogram
   AliRsnListOutput *outMonitorPionDCAz = new AliRsnListOutput("Pion_DCAz", AliRsnListOutput::kHistoDefault);
   outMonitorPionDCAz->AddValue(axisPionDCAz);

   // add outputs to loop
   if (mon) mon->Add(outMonitorPionDCAz);
   if (pdcaz) pdcaz->AddOutput(outMonitorPionDCAz);
  
}

void AddMonitorOutput_PionPIDCut(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *piPID=0)
{

   // Pion PID Cut
   AliRsnValueDaughter *axisPionPIDCut = new AliRsnValueDaughter("pionPID", AliRsnValueDaughter::kTPCnsigmaPi);
   axisPionPIDCut->SetBins(0.0,5,0.01);

   // output: 2D histogram
   AliRsnListOutput *outMonitorPionPIDCut = new AliRsnListOutput("Pion_PID_Cut", AliRsnListOutput::kHistoDefault);
   outMonitorPionPIDCut->AddValue(axisPionPIDCut);

   // add outputs to loop
   if (mon) mon->Add(outMonitorPionPIDCut);
   if (piPID) piPID->AddOutput(outMonitorPionPIDCut);
  
}

void AddMonitorOutput_PionNTPC(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *piNTPC=0)
{

   // Pion PID Cut
   AliRsnValueDaughter *axisPionNTPC = new AliRsnValueDaughter("pionNTPC", AliRsnValueDaughter::kNTPCclusters);
   axisPionNTPC->SetBins(0.0,200,1);

   // output: 2D histogram
   AliRsnListOutput *outMonitorPionNTPC = new AliRsnListOutput("Pion_NTPC", AliRsnListOutput::kHistoDefault);
   outMonitorPionNTPC->AddValue(axisPionNTPC);

   // add outputs to loop
   if (mon) mon->Add(outMonitorPionNTPC);
   if (piNTPC) pNTPC->AddOutput(outMonitorPionNTPC);
  
}

void AddMonitorOutput_PionTPCchi2(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *piTPCchi2=0)
{

   // Pion PID Cut
   AliRsnValueDaughter *axisPionTPCchi2 = new AliRsnValueDaughter("pionTPCchi2", AliRsnValueDaughter::kTPCchi2);
   axisPionTPCchi2->SetBins(0.0,6,.1);

   // output: 2D histogram
   AliRsnListOutput *outMonitorPionTPCchi2 = new AliRsnListOutput("Pion_TPCchi2", AliRsnListOutput::kHistoDefault);
   outMonitorPionTPCchi2->AddValue(axisPionTPCchi2);

   // add outputs to loop
   if (mon) mon->Add(outMonitorPionTPCchi2);
   if (piTPCchi2) pTPCchi2->AddOutput(outMonitorPionTPCchi2);
  
}


void AddMonitorOutput_LambdaP(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lp=0)
{

   // Mass
   AliRsnValueDaughter *axisLambdaP = new AliRsnValueDaughter("lambda_momentum", AliRsnValueDaughter::kP);
   axisLambdaP->SetBins(0.,15.,0.001);

   // output: 2D histogram
   AliRsnListOutput *outMonitorMom = new AliRsnListOutput("Lambda_Momentum", AliRsnListOutput::kHistoDefault);
   outMonitorMom->AddValue(axisLambdaP);

   // add outputs to loop
   if (mon) mon->Add(outMonitorMom);
   if (lp) lp->AddOutput(outMonitorMom);
  
}

void AddMonitorOutput_LambdaPt(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lpt=0)
{

   // Mass
   AliRsnValueDaughter *axisLambdaPt = new AliRsnValueDaughter("lambda_transversemomentum", AliRsnValueDaughter::kV0Pt);
   axisLambdaPt->SetBins(0.,15.,0.001);

   // output: 2D histogram
   AliRsnListOutput *outMonitorTrMom = new AliRsnListOutput("Lambda_TransverseMomentum", AliRsnListOutput::kHistoDefault);
   outMonitorTrMom->AddValue(axisLambdaPt);

   // add outputs to loop
   if (mon) mon->Add(outMonitorTrMom);
   if (lpt) lpt->AddOutput(outMonitorTrMom);
  
}

void AddMonitorOutput_LambdaNegDaughPt(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lnpt=0)
{

   // Mass
   AliRsnValueDaughter *axisLambdaNegDaughPt = new AliRsnValueDaughter("lambda_negdaugh_transversemomentum", AliRsnValueDaughter::kV0NPt);
   axisLambdaNegDaughPt->SetBins(0.,15.,0.001);

   // output: 2D histogram
   AliRsnListOutput *outMonitorLambdaNegDaughTrMom = new AliRsnListOutput("Lambda_NegDaugh_TransverseMomentum", AliRsnListOutput::kHistoDefault);
   outMonitorLambdaNegDaughTrMom->AddValue(axisLambdaNegDaughPt);

   // add outputs to loop
   if (mon) mon->Add(outMonitorLambdaNegDaughTrMom);
   if (lnpt) lnpt->AddOutput(outMonitorLambdaNegDaughTrMom);
  
}

void AddMonitorOutput_LambdaPosDaughPt(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lppt=0)
{

   // Mass
   AliRsnValueDaughter *axisLambdaPosDaughPt = new AliRsnValueDaughter("lambda_posdaugh_transversemomentum", AliRsnValueDaughter::kV0PPt);
   axisLambdaPosDaughPt->SetBins(0.,15.,0.001);

   // output: 2D histogram
   AliRsnListOutput *outMonitorLambdaPosDaughTrMom = new AliRsnListOutput("Lambda_PosDaugh_TransverseMomentum", AliRsnListOutput::kHistoDefault);
   outMonitorLambdaPosDaughTrMom->AddValue(axisLambdaPosDaughPt);

   // add outputs to loop
   if (mon) mon->Add(outMonitorLambdaPosDaughTrMom);
   if (lppt) lppt->AddOutput(outMonitorLambdaPosDaughTrMom);
  
}


void AddMonitorOutput_LambdaMass(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lm=0)
{

   // Mass
   AliRsnValueDaughter *axisMass = new AliRsnValueDaughter("lambda_mass", AliRsnValueDaughter::kV0Mass);
   axisMass->SetBins(1.08,1.16,0.001);

   // output: 2D histogram
   AliRsnListOutput *outMonitorM = new AliRsnListOutput("Lambda_Mass", AliRsnListOutput::kHistoDefault);
   outMonitorM->AddValue(axisMass);

   // add outputs to loop
   if (mon) mon->Add(outMonitorM);
   if (lm) lm->AddOutput(outMonitorM);
  
}

void AddMonitorOutput_LambdaDCA(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *ldca=0)
{
  // Lambda DCA
  AliRsnValueDaughter *axisLambdaDCA = new AliRsnValueDaughter("lambda_dca", AliRsnValueDaughter::kV0DCA);
  axisLambdaDCA->SetBins(0.0,0.4,0.001);
  // output: 2D histogram
  AliRsnListOutput *outMonitorLambdaDCA = new AliRsnListOutput("Lambda_DCA", AliRsnListOutput::kHistoDefault);
  outMonitorLambdaDCA->AddValue(axisLambdaDCA); 
  // add outputs to loop
  if (mon) mon->Add(outMonitorLambdaDCA);
  if (ldca) ldca->AddOutput(outMonitorLambdaDCA);
}

void AddMonitorOutput_LambdaRadius(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *ldca=0)
{
  // Lambda Radius
  AliRsnValueDaughter *axisLambdaRadius = new AliRsnValueDaughter("lambda_radius", AliRsnValueDaughter::kV0Radius);
  axisLambdaRadius->SetBins(0.0,200,0.2);
  // output: 2D histogram
  AliRsnListOutput *outMonitorLambdaRadius = new AliRsnListOutput("Lambda_Radius", AliRsnListOutput::kHistoDefault);
  outMonitorLambdaRadius->AddValue(axisLambdaRadius); 
  // add outputs to loop
  if (mon) mon->Add(outMonitorLambdaRadius);
  if (ldca) ldca->AddOutput(outMonitorLambdaRadius);
}

void AddMonitorOutput_LambdaDaughterDCA(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *ldaugdca=0)
{

   // Lambda Daughter DCA
   AliRsnValueDaughter *axisLambdaDDCA = new AliRsnValueDaughter("lambda_daughterDCA", AliRsnValueDaughter::kDaughterDCA);
   axisLambdaDDCA->SetBins(0.0,2,0.001);

   // output: 2D histogram
   AliRsnListOutput *outMonitorLambdaDDCA = new AliRsnListOutput("Lambda_DaughterDCA", AliRsnListOutput::kHistoDefault);
   outMonitorLambdaDDCA->AddValue(axisLambdaDDCA);

   // add outputs to loop
   if (mon) mon->Add(outMonitorLambdaDDCA);
   if (ldaugdca) ldaugdca->AddOutput(outMonitorLambdaDDCA);
  
}

void AddMonitorOutput_LambdaCosPointAngle(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lcpa=0)
{

   // Lambda Cosine of the Pointing Angle
   AliRsnValueDaughter *axisLambdaCPA = new AliRsnValueDaughter("lambda_cospointang", AliRsnValueDaughter::kCosPointAng);
   axisLambdaCPA->SetBins(0.94,1.,0.0001);

   // output: 2D histogram
   AliRsnListOutput *outMonitorLambdaCPA = new AliRsnListOutput("Lambda_CosineOfPointingAngle", AliRsnListOutput::kHistoDefault);
   outMonitorLambdaCPA->AddValue(axisLambdaCPA);

   // add outputs to loop
   if (mon) mon->Add(outMonitorLambdaCPA);
   if (lcpa) lcpa->AddOutput(outMonitorLambdaCPA);
  
}

void AddMonitorOutput_LambdaProtonPID(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lpPID=0)
{

   // Lambda Cosine of the Pointing Angle
   AliRsnValueDaughter *axisLambdaProtonPID = new AliRsnValueDaughter("lambda_protonPID", AliRsnValueDaughter::kLambdaProtonPIDCut);
   axisLambdaProtonPID->SetBins(0.0,5,0.01);

   // output: 2D histogram
   AliRsnListOutput *outMonitorLambdaProtonPID = new AliRsnListOutput("Lambda_ProtonPID", AliRsnListOutput::kHistoDefault);
   outMonitorLambdaProtonPID->AddValue(axisLambdaProtonPID);

   // add outputs to loop
   if (mon) mon->Add(outMonitorLambdaProtonPID);
   if (lpPID) lpPID->AddOutput(outMonitorLambdaProtonPID);
  
}

void AddMonitorOutput_LambdaPionPID(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lpiPID=0)
{

   // Lambda Cosine of the Pointing Angle
   AliRsnValueDaughter *axisLambdaPionPID = new AliRsnValueDaughter("lambda_pionPID", AliRsnValueDaughter::kLambdaPionPIDCut);
   axisLambdaPionPID->SetBins(0.0,5,0.01);

   // output: 2D histogram
   AliRsnListOutput *outMonitorLambdaPionPID = new AliRsnListOutput("Lambda_PionPID", AliRsnListOutput::kHistoDefault);
   outMonitorLambdaPionPID->AddValue(axisLambdaPionPID);

   // add outputs to loop
   if (mon) mon->Add(outMonitorLambdaPionPID);
   if (lpiPID) lpiPID->AddOutput(outMonitorLambdaPionPID);
  
}

void AddMonitorOutput_LambdaAntiProtonPID(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lapPID=0)
{

   // Lambda Cosine of the Pointing Angle
   AliRsnValueDaughter *axisLambdaAntiProtonPID = new AliRsnValueDaughter("lambda_antiprotonPID", AliRsnValueDaughter::kAntiLambdaAntiProtonPIDCut);
   axisLambdaAntiProtonPID->SetBins(0.0,5,0.01);

   // output: 2D histogram
   AliRsnListOutput *outMonitorLambdaAntiProtonPID = new AliRsnListOutput("Lambda_AntiProtonPID", AliRsnListOutput::kHistoDefault);
   outMonitorLambdaAntiProtonPID->AddValue(axisLambdaAntiProtonPID);

   // add outputs to loop
   if (mon) mon->Add(outMonitorLambdaAntiProtonPID);
   if (lapPID) lapPID->AddOutput(outMonitorLambdaAntiProtonPID);
  
}

void AddMonitorOutput_LambdaAntiPionPID(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lapiPID=0)
{

   // Lambda Cosine of the Pointing Angle
   AliRsnValueDaughter *axisLambdaAntiPionPID = new AliRsnValueDaughter("lambda_antipionPID", AliRsnValueDaughter::kAntiLambdaAntiPionPIDCut);
   axisLambdaAntiPionPID->SetBins(0.0,5,0.01);

   // output: 2D histogram
   AliRsnListOutput *outMonitorLambdaAntiPionPID = new AliRsnListOutput("Lambda_AntiPionPID", AliRsnListOutput::kHistoDefault);
   outMonitorLambdaAntiPionPID->AddValue(axisLambdaAntiPionPID);

   // add outputs to loop
   if (mon) mon->Add(outMonitorLambdaAntiPionPID);
   if (lapiPID) lpiPID->AddOutput(outMonitorLambdaAntiPionPID);
  
}
