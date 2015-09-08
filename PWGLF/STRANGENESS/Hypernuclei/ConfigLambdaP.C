//
// *** Configuration script for LcN*->Lambda-P analysis ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
//
Bool_t ConfigLambdaP
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
   Float_t                 MinDCAToVertexXYlambdaDaugh,
   const char             *suffix,
   AliRsnCutSet           *cutsPair
)
{
   // manage suffix
   if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
   
   // 
   // -- Define track cuts -------------------------------------------------------------------------
   //

   //TString s = ""; s+=trackDCAcut; s+="*(0.0026+0.0050/pt^1.01)";
   //const char *formula = s;
   
   // integrated proton cut --> It should work out of the box
   //  AliRsnCutDaughterSigmaStar2010PP *cutPi = new AliRsnCutDaughterSigmaStar2010PP("cutPionForSigmaStar", AliPID::kProton); //ONLY TPC
   AliRsnCutDaughterLStar2010 *cutPr = new AliRsnCutDaughterLStar2010("cutProtonFromLcN", AliPID::kProton);                 //Should be TPC+TOF
   //cutPr->SetPIDCut(piPIDCut);    // This is effective! fPIDCut used in IsSelected() after the call to cutQuality
   //cutPr->SetMinTPCcluster(NTPCcluster);   // NOTE!!!! Not effective!! fMinTPCclsuter NOT USED in IsSelected()
   AliRsnCutTrackQuality *cutQuality = (AliRsnCutTrackQuality*) cutPr->CutQuality();
   cutQuality->SetAODTestFilterBit(aodFilterBit);
   cutQuality->SetDefaults2011();
   //cutQuality->SetDCARPtFormula(formula);    
   //cutQuality->SetDCARmax(trackDCAcut);	         
   //cutQuality->SetDCARmin(0.01);	         
   

   // cut set
   AliRsnCutSet *cutSetPr = new AliRsnCutSet("setProtonFromLnC", AliRsnTarget::kDaughter);
   cutSetPr->AddCut(cutPr);
   cutSetPr->SetCutScheme(cutPr->GetName());
   // add to task
   Int_t iCutPi = task->AddTrackCuts(cutSetPr);
   
   // quality cuts
   AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("qualityDaughterLambda");   
   esdTrackCuts->SetAcceptKinkDaughters(0); // 0 = kFalse
   esdTrackCuts->SetMaxChi2PerClusterTPC(4);
   esdTrackCuts->SetMinNClustersTPC(NTPCcluster);
   esdTrackCuts->SetRequireTPCRefit();
   //   esdTrackCuts->SetMinDCAToVertexXY(0.05); //-->Ero questo
   //   esdTrackCuts->SetMinDCAToVertexXY(0.15);
   esdTrackCuts->SetMinDCAToVertexXY(MinDCAToVertexXYlambdaDaugh);
   
   // cut lambda
   AliRsnCutV0 *cutLambda = new AliRsnCutV0("cutLambda", kLambda0, AliPID::kProton, AliPID::kPion);
   cutLambda->SetESDtrackCuts(esdTrackCuts);
   cutLambda->SetTolerance(massTol);
   cutLambda->SetMaxDCAVertex(lambdaDCA);
   cutLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
   cutLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
   cutLambda->SetMinTPCcluster(NTPCcluster);
   cutLambda->SetMaxRapidity(0.8);
   // cutLambda->SetAODTestFilterBit(aodFilterBit);
   cutLambda->SetPIDCutProton(pPIDCut);
   cutLambda->SetPIDCutPion(piPIDCut);
   
   // cut set
   AliRsnCutSet *cutSetLambda = new AliRsnCutSet("setLambda", AliRsnTarget::kDaughter);
   cutSetLambda->AddCut(cutLambda);
   cutSetLambda->SetCutScheme(cutLambda->GetName());
   
   // add to task
   Int_t iCutLambda = task->AddTrackCuts(cutSetLambda);
   
   // cut anti-AntiLambda
   AliRsnCutV0 *cutAntiLambda = new AliRsnCutV0("cutAntiLambda", kLambda0Bar, AliPID::kProton, AliPID::kPion);
   cutAntiLambda->SetESDtrackCuts(esdTrackCuts);
   cutAntiLambda->SetTolerance(massTol);
   cutAntiLambda->SetMaxDCAVertex(lambdaDCA);
   cutAntiLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
   cutAntiLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
   cutAntiLambda->SetMinTPCcluster(NTPCcluster);
   cutAntiLambda->SetMaxRapidity(0.8);
   //cutAntiLambda->SetAODTestFilterBit(aodFilterBit);
   cutAntiLambda->SetPIDCutProton(pPIDCut);
   cutAntiLambda->SetPIDCutPion(piPIDCut);
   
   // cut set
   AliRsnCutSet *cutSetAntiLambda = new AliRsnCutSet("setAntiLambda", AliRsnTarget::kDaughter);
   cutSetAntiLambda->AddCut(cutAntiLambda);
   cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
   // add to task
   Int_t iCutAntiLambda = task->AddTrackCuts(cutSetAntiLambda); 
   
   
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

   const Int_t dim = 12;

   Bool_t   use     [dim] = { 1         ,  1            ,  1             ,  1             , 1         ,  1            ,  1             ,  1             , 1         ,  1            ,  1             ,  1             };
   Bool_t   useIM   [dim] = { 1         ,  1            ,  1             ,  1             , 1         ,  1            ,  1             ,  1             , 1         ,  1            ,  1             ,  1             };
   TString  name    [dim] = {"LaPr"     , "ALaAPr"      , "LaAPr"        , "ALaPr"        ,"LaPrMix"  , "ALaAPrMix"   , "LaAPrMix"     , "ALaPrmix"     ,"LaPrRot"  , "ALaAPrRot"   , "LaAPrRot"     , "ALaPrRot"     };
   TString  comp    [dim] = {"PAIR"     , "PAIR"        , "PAIR"         , "PAIR"         ,"MIX"      , "MIX"         , "MIX"          , "MIX"          ,"ROTATE1"  , "ROTATE1"     , "ROTATE1"      , "ROTATE1"      };
   TString  output  [dim] = {"HIST"     , "HIST"        , "HIST"         , "HIST"         ,"HIST"     , "HIST"        , "HIST"         , "HIST"         ,"HIST"     , "HIST"        , "HIST"         , "HIST"         };
   Char_t   charge1 [dim] = {'0'        , '0'           , '0'            , '0'            ,'0'        , '0'           , '0'            , '0'            ,'0'        , '0'           , '0'            , '0'            };
   Char_t   charge2 [dim] = {'+'        , '-'           , '-'            , '+'            ,'+'        , '-'           , '-'            , '+'            ,'+'        , '-'           , '-'            , '+'            };
   Int_t    cutID1  [dim] = {iCutLambda , iCutAntiLambda, iCutLambda     , iCutAntiLambda ,iCutLambda , iCutAntiLambda, iCutLambda     , iCutAntiLambda ,iCutLambda , iCutAntiLambda, iCutLambda     , iCutAntiLambda };
   Int_t    cutID2  [dim] = { iCutPi    ,  iCutPi       , iCutPi         , iCutPi         , iCutPi    ,  iCutPi       , iCutPi         , iCutPi         , iCutPi    ,  iCutPi       , iCutPi         , iCutPi         };
   Int_t    ipdg    [dim] = { 4122112   ,  4122112      , 4122112        , 4122112        , 4122112   ,  4122112      , 4122112        , 4122112        , 4122112   ,  4122112      , 4122112        , 4122112        };
   Double_t mass    [dim] = { 3.235     ,  3.235        , 3.235          , 3.235          , 3.235     ,  3.235        , 3.235          , 3.235          , 3.235     ,  3.235        , 3.235          , 3.235          };
   
   for (Int_t i = 0; i < dim; i++) {
      if (!use[i]) continue;
      if (collSyst) output[i] = "SPARSE";
      // create output
      AliRsnMiniOutput *out = task->CreateOutput(Form("lambdap_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      out->SetCutID(0, cutID1[i]);
      out->SetCutID(1, cutID2[i]);
      out->SetDaughter(0, AliRsnDaughter::kLambda);
      out->SetDaughter(1, AliRsnDaughter::kProton);
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetMotherPDG(ipdg[i]);
      out->SetMotherMass(mass[i]);
      // pair cuts
      out->SetPairCuts(cutsPair);
      // axis X: invmass
      if (useIM[i]) 
         out->AddAxis(imID, 2000, 2, 10.0);
      // axis Y: transverse momentum
      out->AddAxis(ptID, 20, 0.0, 10.0);
      //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
	 
      if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
      
    } 
    
   AddMonitorOutput_ProtonDCA(cutSetPr->GetMonitorOutput());
   AddMonitorOutput_ProtonPIDCut(cutSetPr->GetMonitorOutput());
   AddMonitorOutput_ProtonNTPC(cutSetPr->GetMonitorOutput());
   AddMonitorOutput_ProtonPIDTPC(cutSetPr->GetMonitorOutput());
   AddMonitorOutput_ProtonPIDTOF(cutSetPr->GetMonitorOutput());

   AddMonitorOutput_LambdaMass(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaP(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaPt(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaDCA(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaRadius(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaDaughterDCA(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaCosPointAngle(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaProtonPID(cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaPionPID(cutSetLambda->GetMonitorOutput());

   AddMonitorOutput_LambdaMass(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaP(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaPt(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaDCA(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaRadius(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaDaughterDCA(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaCosPointAngle(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaAntiProtonPID(cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaAntiPionPID(cutSetAntiLambda->GetMonitorOutput());
  

   if (isMC) {
     
     TString mode = "HIST";
     if (collSyst) mode = "SPARSE";
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarP_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kProton);
     out->SetMotherPDG(4122112);
     out->SetMotherMass(3.235);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 3.20, 3.3);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarM_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kProton);
     out->SetMotherPDG(-4122112);
     out->SetMotherMass(3.235);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 3.2, 3.3);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     

     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     
   }
   
   return kTRUE;
}

void AddMonitorOutput_ProtonDCA(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *pdca=0)
{

   // PionDCA
   AliRsnValueDaughter *axisPionDCA = new AliRsnValueDaughter("pion_dca", AliRsnValueDaughter::kDCAXY);
   axisPionDCA->SetBins(-0.5,0.5,0.001);

   // output: 2D histogram
   AliRsnListOutput *outMonitorPionDCA = new AliRsnListOutput("Pion_DCA", AliRsnListOutput::kHistoDefault);
   outMonitorPionDCA->AddValue(axisPionDCA);

   // add outputs to loop
   if (mon) mon->Add(outMonitorPionDCA);
   if (pdca) pdca->AddOutput(outMonitorPionDCA);
  
}

void AddMonitorOutput_ProtonPIDCut(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *piPID=0)
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

void AddMonitorOutput_ProtonNTPC(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *piNTPC=0)
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


void AddMonitorOutput_LambdaP(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lm=0)
{

   // Mass
   AliRsnValueDaughter *axisMass = new AliRsnValueDaughter("lambda_momentum", AliRsnValueDaughter::kP);
   axisMass->SetBins(0.,15.,0.001);

   // output: 2D histogram
   AliRsnListOutput *outMonitorMom = new AliRsnListOutput("Lambda_Momentum", AliRsnListOutput::kHistoDefault);
   outMonitorMom->AddValue(axisMass);

   // add outputs to loop
   if (mon) mon->Add(outMonitorMom);
   if (lm) lm->AddOutput(outMonitorMom);
  
}

void AddMonitorOutput_LambdaPt(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lm=0)
{

   // Mass
   AliRsnValueDaughter *axisMass = new AliRsnValueDaughter("lambda_transversemomentum", AliRsnValueDaughter::kPt);
   axisMass->SetBins(0.,15.,0.001);

   // output: 2D histogram
   AliRsnListOutput *outMonitorTrMom = new AliRsnListOutput("Lambda_TransverseMomentum", AliRsnListOutput::kHistoDefault);
   outMonitorTrMom->AddValue(axisMass);

   // add outputs to loop
   if (mon) mon->Add(outMonitorTrMom);
   if (lm) lm->AddOutput(outMonitorTrMom);
  
}

void AddMonitorOutput_LambdaMass(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lm=0)
{

   // Mass
   AliRsnValueDaughter *axisMass = new AliRsnValueDaughter("lambda_mass", AliRsnValueDaughter::kV0Mass);
   axisMass->SetBins(0.7,1.5,0.001);

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
   axisLambdaCPA->SetBins(0.9,1.,0.0001);

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
//PID
//TPC
void AddMonitorOutput_ProtonPIDTPC(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *pdca=0)
{

  // Momentum
  AliRsnValueDaughter *axisMomTPC = new AliRsnValueDaughter("pTPC", AliRsnValueDaughter::kPtpc);
  axisMomTPC->SetBins(0.0, 10.0, 0.02);

  AliRsnValueDaughter *axisSigTPC = new AliRsnValueDaughter("ProtonTPC", AliRsnValueDaughter::kTPCsignal);
  axisSigTPC->SetBins(0.0, 500.0, 2.0);
 
  AliRsnListOutput *outMonitordEdxTPC = new AliRsnListOutput("dEdx_VsPtpc_Proton", AliRsnListOutput::kHistoDefault);
  outMonitordEdxTPC->AddValue(axisMomTPC);
  outMonitordEdxTPC->AddValue(axisSigTPC);
  if (mon) mon->Add(outMonitordEdxTPC);
  if (pdca) pdca->AddOutput(outMonitordEdxTPC);

}

//TOF

void AddMonitorOutput_ProtonPIDTOF(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lm=0){

  AliRsnValueDaughter *axisTOFnsigmaP = new AliRsnValueDaughter("p1", AliRsnValueDaughter::kTOFnsigmaP);
  axisTOFnsigmaP->SetBins(-10.,10., 0.1);

  // AliRsnValueDaughter *axisTPCnsigmaP = new AliRsnValueDaughter("p2", AliRsnValueDaughter::kTPCnsigmaP);
  // axisTPCnsigmaP->SetBins(-10.,10., 0.1);

  // kTOFdeltaP
  AliRsnValueDaughter *axisTOFdeltaP = new AliRsnValueDaughter("Dp", AliRsnValueDaughter::kTOFdeltaP);
  axisTOFdeltaP->SetBins(-3000.,3000., 10.);

  AliRsnValueDaughter *axisMomP = new AliRsnValueDaughter("p3", AliRsnValueDaughter::kP);
  axisMomP->SetBins(0.0, 10.0, 0.02);

  AliRsnListOutput *outMonitorTOFnsigmaP = new AliRsnListOutput("TOF_nsigmaPro_vsP", AliRsnListOutput::kHistoDefault);
  outMonitorTOFnsigmaP->AddValue(axisMomP);
  outMonitorTOFnsigmaP->AddValue(axisTOFnsigmaP);
  if (mon) mon->Add(outMonitorTOFnsigmaP);
  if (lm) lm->AddOutput(outMonitorTOFnsigmaP);

  // output: 2D histogram of TOF signal vs. TOF momentum
  AliRsnListOutput *outMonitorTOFdeltaP = new AliRsnListOutput("TOF_deltaPro_vsP", AliRsnListOutput::kHistoDefault);
  outMonitorTOFdeltaP->AddValue(axisMomP);
  outMonitorTOFdeltaP->AddValue(axisTOFdeltaP);
  if (mon) mon->Add(outMonitorTOFdeltaP);
  if (lm) lm->AddOutput(outMonitorTOFdeltaP);

}
