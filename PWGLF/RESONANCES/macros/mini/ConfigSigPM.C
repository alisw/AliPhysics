/***************************************************************************
Paraskevi.Ganoti@cern.ch - created on 2019 with the help of F. Bellini and 
other scripts in the repo.  
Configuration script for SigmaStar(9801385) analysis
****************************************************************************/

Bool_t ConfigSigPM(AliRsnMiniAnalysisTask *task, 
		Bool_t                 isMC, 
		AliPIDResponse::EBeamType collSys = AliPIDResponse::kPBPB, //=0, kPPB=1, kPBPB=2
		AliRsnCutSet           *cutsPair,             //cuts on the pair
                AliRsnCutSet           *cutsPairY,             //cuts on the pair Y only
		Bool_t                 enaMultSel = kTRUE,    //enable multiplicity axis
      		Float_t                masslow = 1.2,         //inv mass axis low edge 
		Float_t                massup = 2.,          //inv mass axis upper edge 
		Int_t                  nbins = 800,           //inv mass axis n bins
		Float_t                nsigma = 3.0,          //nsigma of TPC PID cut
		Bool_t                 enableMonitor = kTRUE, //enable single track QA plots
                Float_t                pi_Ls_PIDCut=4.,        //nsigma V0 daughters
                Float_t                LsDCA = 0.3,             //V0 vtx to PV
                Float_t                LsCosPoinAn = 0.98,      // cos of Pointing Angle
                Float_t                LsDaughDCA=0.8,          // dca V0 daughters
                Float_t                massTol = 0.006,         // mass tolerance 6 MeV
                Float_t                massTolVeto = 0.0043,    // mass tol veto
                Bool_t                 Switch = kFALSE,         // switch
                Float_t                pLife = 25,              // life
                Float_t                v0rapidity= 0.5,         // V0 rapidity
                Float_t                radiuslow=5.,            // radius low
                Bool_t                 doCustomDCAcuts=kTRUE,   //custom dca cuts for V0 daughters
                Double_t               dcaProton=0.1,           // proton dca
		Double_t               dcaPion=0.1,             //pion dca
		Int_t                  pidCUT=1)             //pion PID cut set, 1 for nominal, 2 for systematic check
{
  //-----------------------
  //General 
  //-----------------------
  TString partname="SigmaStar";
//  RSNPID  d1 = AliRsnDaughter::kPion;
//  RSNPID  d2 = AliRsnDaughter::kLambda;
  Int_t   aodFilterBit = 0;

  //Additional options for monitoring plots
  TString monitorOpt = "NoSIGN";
  
  //-----------------------
  // CUTS
  //-----------------------
  //use default quality cuts std 2010 with crossed rows TPC
  Bool_t useCrossedRows = 1;

  AliRsnCutTrackQuality * fCutQuality = new AliRsnCutTrackQuality("CutQuality");
  fCutQuality->SetDefaults2011(useCrossedRows, kFALSE);
  
  AliRsnCutTOFMatch  *iCutTOFMatch     = new AliRsnCutTOFMatch("CutTOFMatch");
  AliRsnCutPIDNSigma *iCutTPCNSigma    = new AliRsnCutPIDNSigma("CutTPCNSigma", AliPID::kPion, AliRsnCutPIDNSigma::kTPC);//, AliRsnCutPIDNSigma::kTPCinnerP );
  AliRsnCutPIDNSigma *iCutTOFNSigma    = new AliRsnCutPIDNSigma("CutTOFNSigma", AliPID::kPion, AliRsnCutPIDNSigma::kTOF);//, AliRsnCutPIDNSigma::kP );
  AliRsnCutPIDNSigma *iCutTPCTOFNSigma = new AliRsnCutPIDNSigma("CutTPCTOFNSigma", AliPID::kPion, AliRsnCutPIDNSigma::kTPC);
  //for setting PID cuts without any selection in pT
  //iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);

  if (pidCUT==1) {
    //for setting PID cuts in given pT ranges
    iCutTPCNSigma->AddPIDRange(5.0, 0.0, 0.35);
    iCutTPCNSigma->AddPIDRange(3.0, 0.35, 0.5);
    iCutTPCNSigma->AddPIDRange(2.0, 0.5, 20.);
    
    iCutTOFNSigma->AddPIDRange(3.0, 0.0, 1.5);
    iCutTOFNSigma->AddPIDRange(2.5, 1.5, 1E20);

     AliRsnCutSet * myCutSet = new AliRsnCutSet("MyCutSet", AliRsnTarget::kDaughter);
     myCutSet->AddCut(fCutQuality);
     myCutSet->AddCut(iCutTPCNSigma);
     myCutSet->AddCut(iCutTOFMatch);
     myCutSet->AddCut(iCutTOFNSigma);
     // myCutSet->AddCut(iCutTPCTOFNSigma); 
     
     myCutSet->SetCutScheme( Form("%s & ((%s & (!%s))|(%s&%s))",fCutQuality->GetName(), iCutTPCNSigma->GetName(), iCutTOFMatch->GetName(),iCutTOFNSigma->GetName(), iCutTPCNSigma->GetName()) ) ;
    
  }

  if (pidCUT==2) {
    //for setting PID cuts in given pT ranges
    iCutTPCNSigma->AddPIDRange(3.0, 0.0, 0.5);
    iCutTPCNSigma->AddPIDRange(2.0, 0.5, 1E20);
    iCutTPCTOFNSigma->SinglePIDRange(5.0);
    iCutTOFNSigma->AddPIDRange(3.0, 0.0, 1E20);

     AliRsnCutSet * myCutSet = new AliRsnCutSet("MyCutSet", AliRsnTarget::kDaughter);
     myCutSet->AddCut(fCutQuality);
     myCutSet->AddCut(iCutTPCNSigma);
     myCutSet->AddCut(iCutTOFMatch);
     myCutSet->AddCut(iCutTOFNSigma);
     myCutSet->AddCut(iCutTPCTOFNSigma);

     // scheme:
     // quality & [ (TOF & TPCTOF) || (!TOFmatch & TPConly) ]
     myCutSet->SetCutScheme(Form("%s&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName())) ;
  }

  Int_t icutPi = task->AddTrackCuts(myCutSet);
 
  //set daughter cuts
  Int_t icut1 = icutPi;
  //Int_t icut2 = icutPi;

   // selections for daugthers of Lambdas

  Float_t crossedRows = 70;
  Float_t rowsbycluster = 0.8;
  Float_t DCAxy = 0.1; //0.06;
     
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("qualityDaughterLs");
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.15,1.E10);
  esdTrackCuts->SetRequireTPCRefit();
  esdTrackCuts->SetAcceptKinkDaughters(0); //
  //esdTrackCuts->SetMinNCrossedRowsTPC(useCrossedRows);
  //esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(rowsbycluster);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetMinNClustersTPC(70);
//  esdTrackCuts->SetMinDCAToVertexXY(DCAxy); //Use one of the two - pt dependent or fixed value cut.
  
  //V0s

  // selections for Lambdas
    AliRsnCutV0 *cutLambdas = new AliRsnCutV0("cutLambdas", kLambda0, AliPID::kProton, AliPID::kPion);
    cutLambdas->SetPIDCutPion(pi_Ls_PIDCut);        // PID for the pion daughter of Ls
    cutLambdas->SetPIDCutProton(pi_Ls_PIDCut);        // PID for the proton daughter of Ls
    cutLambdas->SetESDtrackCuts(esdTrackCuts);
    cutLambdas->SetMaxDaughtersDCA(LsDaughDCA);
    //cutLambdas->SetMaxDCAVertex(LsDCA);
    cutLambdas->SetMinCosPointingAngle(LsCosPoinAn);
    cutLambdas->SetTolerance(massTol);
    //cutLambdas->SetToleranceVeto(massTolVeto);   //Rejection range for Competing V0 Rejection
    cutLambdas->SetSwitch(Switch);
    cutLambdas->SetfLife(pLife);
    cutLambdas->SetfLowRadius(radiuslow);
    cutLambdas->SetfHighRadius(200);
    cutLambdas->SetMaxRapidity(v0rapidity);
    cutLambdas->SetDifferentDCACutPosNegTrack(doCustomDCAcuts);
    cutLambdas->SetMinDCAToVtxXYPositiveTrack(dcaProton);
    cutLambdas->SetMinDCAToVtxXYNegativeTrack(dcaPion);
    // cutLambdas->SetpT_Tolerance(tol_switch);
    //cutLambdas->SetMassTolSigma(tol_sigma);

    AliRsnCutSet *cutSetLs = new AliRsnCutSet("setLs", AliRsnTarget::kDaughter);
    cutSetLs->AddCut(cutLambdas);
    cutSetLs->SetCutScheme(cutLambdas->GetName());
    Int_t iCutLs = task->AddTrackCuts(cutSetLs);

    AliRsnCutV0 *cutAntiLambdas = new AliRsnCutV0("cutAntiLambdas", kLambda0Bar, AliPID::kProton, AliPID::kPion);
    cutAntiLambdas->SetPIDCutPion(pi_Ls_PIDCut);        // PID for the pion daughter of Ls
    cutAntiLambdas->SetPIDCutProton(pi_Ls_PIDCut);        // PID for the proton daughter of Ls
    cutAntiLambdas->SetESDtrackCuts(esdTrackCuts);
    cutAntiLambdas->SetMaxDaughtersDCA(LsDaughDCA);
//    cutAntiLambdas->SetMaxDCAVertex(LsDCA);
    cutAntiLambdas->SetMinCosPointingAngle(LsCosPoinAn);
    cutAntiLambdas->SetTolerance(massTol);
    //cutLambdas->SetToleranceVeto(massTolVeto);   //Rejection range for Competing V0 Rejection
    cutAntiLambdas->SetSwitch(Switch);
    cutAntiLambdas->SetfLife(pLife);
    cutAntiLambdas->SetfLowRadius(radiuslow);
    cutAntiLambdas->SetfHighRadius(200);
    cutAntiLambdas->SetMaxRapidity(v0rapidity);
    cutAntiLambdas->SetDifferentDCACutPosNegTrack(doCustomDCAcuts);
    cutAntiLambdas->SetMinDCAToVtxXYPositiveTrack(dcaPion);
    cutAntiLambdas->SetMinDCAToVtxXYNegativeTrack(dcaProton);
    // cutLambdas->SetpT_Tolerance(tol_switch);
    //cutLambdas->SetMassTolSigma(tol_sigma);

    AliRsnCutSet *cutSetAntiLs = new AliRsnCutSet("setAntiLs", AliRsnTarget::kDaughter);
    cutSetAntiLs->AddCut(cutAntiLambdas);
    cutSetAntiLs->SetCutScheme(cutAntiLambdas->GetName());
    Int_t iCutAntiLs = task->AddTrackCuts(cutSetAntiLs);

  // end V0s

  //QA plots 
  if (enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
//    AddMonitorOutput(isMC, cutSetQuality->GetMonitorOutput(), monitorOpt.Data());
      AddMonitorOutput(isMC, myCutSet->GetMonitorOutput(), monitorOpt.Data(), 0);
//    AddMonitorOutput(isMC, cutSetAntiLs->GetMonitorOutput(), monitorOpt.Data(), 0);
  }  

  //-----------------------
  // OUTPUT
  //-----------------------
  /* invariant mass       */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM resolution        */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum     */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality           */ Int_t multID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity       */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity             */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
  /* 1st daughter pt      */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt, kFALSE);
  /* 2nd daughter pt      */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt, kFALSE);
  /* 1st daughter p       */ Int_t fdp    = task->CreateValue(AliRsnMiniValue::kFirstDaughterP, kFALSE);
  /* 2nd daughter p       */ Int_t sdp    = task->CreateValue(AliRsnMiniValue::kSecondDaughterP, kFALSE);
  /* Asymmetry data       */ Int_t asymd   = task->CreateValue(AliRsnMiniValue::kAsym, kFALSE);
  /* Asymmetry            */ Int_t asym   = task->CreateValue(AliRsnMiniValue::kAsym, kTRUE);
   /* gener. transv. mom.  */
  if (isMC) Int_t ptIDgen= task->CreateValue(AliRsnMiniValue::kPt, kTRUE);
  
  // TString output[8] = {"HIST",      "HIST",      "HIST",      "HIST",       "HIST",      "HIST",      "HIST",      "HIST"}; // or "SPARSE"
  // TString name[8] =   {"SigmaP",   "ASigmaM",   "ASigmaP",   "ASigmaM",   "SigmaPmix",  "SigmaMmix", "ASigmaPmix", "ASigmaMmix" };
  // TString comp[8] =   {"PAIR",      "PAIR",       "PAIR",      "PAIR",       "MIX",        "MIX",       "MIX",        "MIX"  };
  // Char_t charge1[8] = {  '0',         '0',          '0',         '0',          '0',         '0',         '0',           '0' };
  // Char_t charge2[8] = {  '+',         '-',          '-',         '+',          '+',         '-',         '-',           '+'};
  // Int_t ipdg [8] =    {  3224,        3114,        -3224,       -3114,         3224,        3114,       -3224,          -3114 };
  // Double_t mass [8]=  { 1.3828,       1.3872,      1.3828,      1.3872,       1.3828,      1.3872,      1.3828,         1.3872}; 

  TString output[8] = {"HIST",      "HIST",          "HIST",      "HIST",          "HIST",        "HIST",        "HIST",           "HIST"    }; // or "SPARSE"
  TString name[8] =   {"SigmaP",   "SigmaM",     "SigmaPmix",    "SigmaMmix",    "ASigmaP",     "ASigmaM",    "ASigmaPmix",     "ASigmaMmix" };
  TString comp[8] =   {"PAIR",      "PAIR",        "MIX",         "MIX",           "PAIR",       "PAIR",         "MIX",             "MIX"    };
  Char_t charge1[8] = {  '0',         '0',          '0',            '0',             '0',          '0',           '0',               '0'     };
  Char_t charge2[8] = {  '+',         '-',          '+',            '-',             '-',          '+',           '-',               '+'     };
  Int_t ipdg [8] =    {  3224,        3114,          3224,          3114,           -3224 ,       -3114,         -3224 ,            -3114    };
  Double_t mass [8]=  { 1.3828,       1.3872,        1.3828,       1.3872,          1.3828,       1.3872,        1.3828,            1.3872   }; 
  Int_t cutID1 [8]= {iCutLs, iCutLs, iCutLs, iCutLs, iCutAntiLs, iCutAntiLs, iCutAntiLs, iCutAntiLs}; 
  Int_t cutID2 [8]= {icutPi, icutPi, icutPi, icutPi, icutPi, icutPi, icutPi, icutPi};

  //DATA 
  for (Int_t i = 0; i < 8; i++) {
    output[i] = "SPARSE";
    AliRsnMiniOutput *out = task->CreateOutput(Form("SigPM_%s", name[i].Data()), output[i].Data(), comp[i].Data());
    out->SetCutID(0, cutID1[i]);
    out->SetCutID(1, cutID2[i]);
    out->SetDaughter(0, AliRsnDaughter::kLambda);
    out->SetDaughter(1, AliRsnDaughter::kPion);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(ipdg[i]);
    out->SetMotherMass(mass[i]);
    out->SetPairCuts(cutsPair);
    // axis X: invmass 
    out->AddAxis(imID, nbins, masslow, massup);
    //axis Y: mother pt
    out->AddAxis(ptID, 150, 0.0, 15.0); //default use mother pt
    //axis Z: multiplicity
      if (enaMultSel) out->AddAxis(multID, 100, 0.0, 100.0);
  }
  
    AliRsnMiniOutput *outasmd = task->CreateOutput("hAsymmetryDAtaSp", "HIST", "PAIR");
    outasmd->SetDaughter(0, AliRsnDaughter::kLambda);
    outasmd->SetDaughter(1, AliRsnDaughter::kPion);
    outasmd->SetCutID(0, iCutLs);
    outasmd->SetCutID(1, icutPi);
    outasmd->SetMotherPDG(3224);  
    outasmd->SetMotherMass(1.3828);
    outasmd->SetPairCuts(cutsPair);
    outasmd->SetDaughter(0, AliRsnDaughter::kLambda);
    outasmd->SetDaughter(1, AliRsnDaughter::kPion);
    outasmd->SetCharge(0, charge1[0]);
    outasmd->SetCharge(1, charge2[0]);
    outasmd->AddAxis(asymd, 200, -1.0, 1.0);
    if (enaMultSel) outasmd->AddAxis(multID, 100, 0.0, 100.0);

   //AliRsnMiniOutput *out = task->CreateOutput("hOpeningAngle", "HIST", "PAIR");
   //out->SetCutID(0, cutID1[1]);
   //out->SetCutID(1, cutID2[1]);
   //out->SetPairCuts(cutsPair);
   //out->AddAxis(ptID, 350, 0.0, 35.0);
   //out->AddAxis(opAngl, 300, 0., 150.);
 
   AddMonitorOutput_LambdaPt(cutSetLs->GetMonitorOutput());
   AddMonitorOutput_LambdaNegDaughPt(cutSetLs->GetMonitorOutput());
   AddMonitorOutput_LambdaPosDaughPt(cutSetLs->GetMonitorOutput());
   AddMonitorOutput_LambdaMass(cutSetLs->GetMonitorOutput());
   AddMonitorOutput_LambdaDCA(cutSetLs->GetMonitorOutput());
   AddMonitorOutput_LambdaRadius(cutSetLs->GetMonitorOutput());
   AddMonitorOutput_LambdaDaughterDCA(cutSetLs->GetMonitorOutput());
   AddMonitorOutput_LambdaCosPointAngle(cutSetLs->GetMonitorOutput());
   AddMonitorOutput_LambdaProtonPID(cutSetLs->GetMonitorOutput());
   AddMonitorOutput_LambdaPionPID(cutSetLs->GetMonitorOutput());

   AddMonitorOutput_LambdaPt(cutSetAntiLs->GetMonitorOutput());
   AddMonitorOutput_LambdaNegDaughPt(cutSetAntiLs->GetMonitorOutput());
   AddMonitorOutput_LambdaPosDaughPt(cutSetAntiLs->GetMonitorOutput());
   AddMonitorOutput_LambdaMass(cutSetAntiLs->GetMonitorOutput());
   AddMonitorOutput_LambdaDCA(cutSetAntiLs->GetMonitorOutput());
   AddMonitorOutput_LambdaRadius(cutSetAntiLs->GetMonitorOutput());
   AddMonitorOutput_LambdaDaughterDCA(cutSetAntiLs->GetMonitorOutput());
   AddMonitorOutput_LambdaCosPointAngle(cutSetAntiLs->GetMonitorOutput());
   AddMonitorOutput_LambdaAntiProtonPID(cutSetAntiLs->GetMonitorOutput());
   AddMonitorOutput_LambdaAntiPionPID(cutSetAntiLs->GetMonitorOutput());
   
   
 //   AddMonitorOutput_LambdaMass(cutSetAntiLambda->GetMonitorOutput());
//    //AddMonitorOutput_LambdaP(cutSetAntiLambda->GetMonitorOutput());
//    AddMonitorOutput_LambdaPt(cutSetAntiLambda->GetMonitorOutput());
//    AddMonitorOutput_LambdaNegDaughPt(cutSetAntiLambda->GetMonitorOutput());
//    AddMonitorOutput_LambdaPosDaughPt(cutSetAntiLambda->GetMonitorOutput());
//    AddMonitorOutput_LambdaDCA(cutSetAntiLambda->GetMonitorOutput());
//    AddMonitorOutput_LambdaRadius(cutSetAntiLambda->GetMonitorOutput());
//    AddMonitorOutput_LambdaDaughterDCA(cutSetAntiLambda->GetMonitorOutput());
//    AddMonitorOutput_LambdaCosPointAngle(cutSetAntiLambda->GetMonitorOutput());
//    AddMonitorOutput_LambdaAntiProtonPID(cutSetAntiLambda->GetMonitorOutput());
// AddMonitorOutput_LambdaAntiPionPID(cutSetAntiLambda->GetMonitorOutput());
  // //Template for BG
  // TString bgTemplate[7]  = {"rho", "omega", "Kstar", "antiKstar", "K0s", "phi","f2"};
  // Char_t bgTemplateC1[7] = {'+', '+', '+', '+', '+', '+', '+'};
  // Char_t bgTemplateC2[7] = {'-', '-', '-', '-', '-', '-', '-'};
  // Int_t bgTemplatePDG[7] = {113, 223, 313, -313, 310, 333, 225};
  // Int_t bgTemplateM[7]   = {775.26, 8.49, 891.66, 891.66, 497.611, 1019.461, 1275.5};
  // RSNPID bgID1[7] = {AliRsnDaughter::kPion, AliRsnDaughter::kPion,AliRsnDaughter::kKaon, AliRsnDaughter::kPion, AliRsnDaughter::kPion, AliRsnDaughter::kKaon, AliRsnDaughter::kPion};
  // RSNPID bgID2[7] = {AliRsnDaughter::kPion, AliRsnDaughter::kPion,AliRsnDaughter::kPion, AliRsnDaughter::kKaon, AliRsnDaughter::kPion, AliRsnDaughter::kKaon, AliRsnDaughter::kPion};

  if (isMC) {

    //TRUE RECO PAIRS - MASS
    output[0] = "SPARSE";
    // create output
    AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarP_TrueMC_%s", name[0].Data()), output[0].Data(), "TRUE");
    // selection settings
    out->SetCutID(0, cutID1[0]);
    out->SetCutID(1, cutID2[0]);
    out->SetDaughter(0, AliRsnDaughter::kLambda);
    out->SetDaughter(1, AliRsnDaughter::kPion);
    out->SetCharge(0, charge1[0]);
    out->SetCharge(1, charge2[0]);
    out->SetMotherPDG(ipdg[0]);
    out->SetMotherMass(mass[0]);
    out->SetPairCuts(cutsPair);
    // binnings
    out->AddAxis(imID, 800, 1.2, 2.0);
    out->AddAxis(ptID, 150, 0.0, 15.0);
    //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
    
    if (enaMultSel) out->AddAxis(multID, 100, 0.0, 100.0);

    output[1] = "SPARSE";
    // create output
    AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarM_TrueMC_%s", name[1].Data()), output[1].Data(), "TRUE");
    // selection settings
    out->SetCutID(0, cutID1[1]);
    out->SetCutID(1, cutID2[1]);
    out->SetDaughter(0, AliRsnDaughter::kLambda);
    out->SetDaughter(1, AliRsnDaughter::kPion);
    out->SetCharge(0, charge1[1]);
    out->SetCharge(1, charge2[1]);
    out->SetMotherPDG(ipdg[1]);
    out->SetMotherMass(mass[1]);
    out->SetPairCuts(cutsPair);
    // binnings
    out->AddAxis(imID, 800, 1.2, 2.0);
    out->AddAxis(ptID, 150, 0.0, 15.0);
    //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
    
    if (enaMultSel) out->AddAxis(multID, 100, 0.0, 100.0);

    output[4] = "SPARSE";
    // create output
    AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarPbar_TrueMC_%s", name[4].Data()), output[4].Data(), "TRUE");
    // selection settings
    out->SetCutID(0, cutID1[4]);
    out->SetCutID(1, cutID2[4]);
    out->SetDaughter(0, AliRsnDaughter::kLambda);
    out->SetDaughter(1, AliRsnDaughter::kPion);
    out->SetCharge(0, charge1[4]);
    out->SetCharge(1, charge2[4]);
    out->SetMotherPDG(ipdg[4]);
    out->SetMotherMass(mass[4]);
    out->SetPairCuts(cutsPair);
    // binnings
    out->AddAxis(imID, 800, 1.2, 2.0);
    out->AddAxis(ptID, 150, 0.0, 15.0);
    //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
    
    if (enaMultSel) out->AddAxis(multID, 100, 0.0, 100.0);

    output[5] = "SPARSE";
    // create output
    AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarMbar_TrueMC_%s", name[5].Data()), output[5].Data(), "TRUE");
    // selection settings
    out->SetCutID(0, cutID1[5]);
    out->SetCutID(1, cutID2[5]);
    out->SetDaughter(0, AliRsnDaughter::kLambda);
    out->SetDaughter(1, AliRsnDaughter::kPion);
    out->SetCharge(0, charge1[5]);
    out->SetCharge(1, charge2[5]);
    out->SetMotherPDG(ipdg[5]);
    out->SetMotherMass(mass[5]);
    out->SetPairCuts(cutsPair);
    // binnings
    out->AddAxis(imID, 800, 1.2, 2.0);
    out->AddAxis(ptID, 150, 0.0, 15.0);

    if (enaMultSel) out->AddAxis(multID, 100, 0.0, 100.0);

   //AliRsnMiniOutput *out = task->CreateOutput("hAsymmetryMC", "HIST", "TRUE");
   // out->SetCutID(0, cutID1[0]);
   // out->SetCutID(1, cutID2[0]);
   // out->SetDaughter(0, AliRsnDaughter::kLambda);
   // out->SetDaughter(1, AliRsnDaughter::kPion);
   // out->SetCharge(0, charge1[0]);
   // out->SetCharge(1, charge2[0]);
   // out->SetMotherPDG(ipdg[0]); 
   // out->SetMotherMass(mass[0]);
   // out->SetPairCuts(cutsPair);
   // out->AddAxis(ptID, 100, 0.0, 10.0);
   // out->AddAxis(opAngl, 300, 0., 150.);

   // AliRsnMiniOutput *out = task->CreateOutput("hOpeningAngleMCTrueM", "HIST", "TRUE");
   // out->SetCutID(0, cutID1[1]);
   // out->SetCutID(1, cutID2[1]);
   // out->SetDaughter(0, AliRsnDaughter::kLambda);
   // out->SetDaughter(1, AliRsnDaughter::kPion);
   // out->SetCharge(0, charge1[1]);
   // out->SetCharge(1, charge2[1]);
   // out->SetMotherPDG(ipdg[1]);
   // out->SetMotherMass(mass[1]);
   // out->SetPairCuts(cutsPair);
   // out->AddAxis(ptID, 100, 0.0, 10.0);
   // out->AddAxis(opAngl, 300, 0., 150.);
 
   // AliRsnMiniOutput *out = task->CreateOutput("hOpeningAngleMCTrueAP", "HIST", "TRUE");
   // out->SetCutID(0, cutID1[4]);
   // out->SetCutID(1, cutID2[4]);
   // out->SetDaughter(0, AliRsnDaughter::kLambda);
   // out->SetDaughter(1, AliRsnDaughter::kPion);
   // out->SetCharge(0, charge1[4]);
   // out->SetCharge(1, charge2[4]);
   // out->SetMotherPDG(ipdg[4]);
   // out->SetMotherMass(mass[4]);
   // out->SetPairCuts(cutsPair);
   // out->AddAxis(ptID, 100, 0.0, 10.0);
   // out->AddAxis(opAngl, 300, 0., 150.);
    
   // AliRsnMiniOutput *out = task->CreateOutput("hOpeningAngleMCTrueAM", "HIST", "TRUE");
   // out->SetCutID(0, cutID1[5]);
   // out->SetCutID(1, cutID2[5]);
   // out->SetDaughter(0, AliRsnDaughter::kLambda);
   // out->SetDaughter(1, AliRsnDaughter::kPion);
   // out->SetCharge(0, charge1[5]);
   // out->SetCharge(1, charge2[5]);
   // out->SetMotherPDG(ipdg[5]);
   // out->SetMotherMass(mass[5]);
   // out->SetPairCuts(cutsPair);
   // out->AddAxis(ptID, 100, 0.0, 10.0);
   // out->AddAxis(opAngl, 300, 0., 150.);
    
    //GENERATED PAIRS
    AliRsnMiniOutput *outasm = task->CreateOutput("hAsymmetryMCSp", "HIST", "MOTHER");
    outasm->SetDaughter(0, AliRsnDaughter::kLambda);
    outasm->SetDaughter(1, AliRsnDaughter::kPion);
    outasm->SetMotherPDG(3224);
    outasm->SetMotherMass(1.3828);
    outasm->SetPairCuts(cutsPairY);
    outasm->AddAxis(asym, 200, -1.0, 1.0);
    if (enaMultSel) outasm->AddAxis(multID, 100, 0.0, 100.0);

    AliRsnMiniOutput *outasm = task->CreateOutput("hAsymmetryMCSm", "HIST", "MOTHER");
    outasm->SetDaughter(0, AliRsnDaughter::kLambda);
    outasm->SetDaughter(1, AliRsnDaughter::kPion);
    outasm->SetMotherPDG(3114);   
    outasm->SetMotherMass(1.3872);
    outasm->SetPairCuts(cutsPairY);
    outasm->AddAxis(asym, 200, -1.0, 1.0);
    if (enaMultSel) outasm->AddAxis(multID, 100, 0.0, 100.0);

    AliRsnMiniOutput *outasm = task->CreateOutput("hAsymmetryMCSpbar", "HIST", "MOTHER");
    outasm->SetDaughter(0, AliRsnDaughter::kLambda);
    outasm->SetDaughter(1, AliRsnDaughter::kPion);
    outasm->SetMotherPDG(-3224);
    outasm->SetMotherMass(1.3828); 
    outasm->SetPairCuts(cutsPairY);
    outasm->AddAxis(asym, 200, -1.0, 1.0);
    if (enaMultSel) outasm->AddAxis(multID, 100, 0.0, 100.0);

    AliRsnMiniOutput *outasm = task->CreateOutput("hAsymmetryMCSmbar", "HIST", "MOTHER");
    outasm->SetDaughter(0, AliRsnDaughter::kLambda);
    outasm->SetDaughter(1, AliRsnDaughter::kPion);
    outasm->SetMotherPDG(-3114);
    outasm->SetMotherMass(1.3872); 
    outasm->SetPairCuts(cutsPairY);
    outasm->AddAxis(asym, 200, -1.0, 1.0);
    if (enaMultSel) outasm->AddAxis(multID, 100, 0.0, 100.0);

    output[0] = "HIST";
    AliRsnMiniOutput * outm = task->CreateOutput(Form("motherSigmaPpt_%s", name[0].Data()), output[0].Data(),"MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kLambda);
    outm->SetDaughter(1, AliRsnDaughter::kPion);
    outm->SetMotherPDG(3224);
    outm->SetMotherMass(1.3828);
    outm->SetPairCuts(cutsPairY);
    outm->AddAxis(ptIDgen, 150, 0.0, 15.0);
    if (enaMultSel) outm->AddAxis(multID, 100, 0.0, 100.0);

    output[1] = "HIST";
    AliRsnMiniOutput * outm = task->CreateOutput(Form("motherSigmaMpt_%s", name[1].Data()), output[1].Data(),"MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kLambda);
    outm->SetDaughter(1, AliRsnDaughter::kPion);
    outm->SetMotherPDG(3114);
    outm->SetMotherMass(1.3872);
    outm->SetPairCuts(cutsPairY);
    outm->AddAxis(ptIDgen, 150, 0.0, 15.0);
    if (enaMultSel) outm->AddAxis(multID, 100, 0.0, 100.0);

    output[4] = "HIST";
    AliRsnMiniOutput * outm = task->CreateOutput(Form("motherSigmaPbarpt_%s", name[4].Data()), output[4].Data(),"MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kLambda);
    outm->SetDaughter(1, AliRsnDaughter::kPion);
    outm->SetMotherPDG(-3224);
    outm->SetMotherMass(1.3828);
    outm->SetPairCuts(cutsPairY);
    outm->AddAxis(ptIDgen, 150, 0.0, 15.0);
    if (enaMultSel) outm->AddAxis(multID, 100, 0.0, 100.0);

    output[5] = "HIST";
    AliRsnMiniOutput * outm = task->CreateOutput(Form("motherSigmaMbarpt_%s", name[5].Data()), output[5].Data(),"MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kLambda);
    outm->SetDaughter(1, AliRsnDaughter::kPion);
    outm->SetMotherPDG(-3114);
    outm->SetMotherMass(1.3872);
    outm->SetPairCuts(cutsPairY);
    outm->AddAxis(ptIDgen, 150, 0.0, 15.0);
    if (enaMultSel) outm->AddAxis(multID, 100, 0.0, 100.0);

//    AliRsnMiniOutput *out = task->CreateOutput("hOpeningAngleGener", "HIST", "PAIR");
//    outm->SetDaughter(0, AliRsnDaughter::kLambda);
//    outm->SetDaughter(1, AliRsnDaughter::kPion);
//    outm->SetMotherPDG(-3114);  
//    outm->SetMotherMass(1.3872);
//    out->SetPairCuts(cutsPair);
//    out->AddAxis(ptIDgen, 100, 0.0, 10.0);
//    out->AddAxis(opAnglsim, 300, 0., 150.);
 
  }

  return kTRUE;
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
   if (lapiPID) lapiPID->AddOutput(outMonitorLambdaAntiPionPID);
  
}

