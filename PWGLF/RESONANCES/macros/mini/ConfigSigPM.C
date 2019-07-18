/***************************************************************************
Paraskevi.Ganoti@cern.ch - created on 2019 with the help of F. Bellini and 
other scripts in the repo.  
Configuration script for SigmaStar(9801385) analysis
****************************************************************************/

Bool_t ConfigSigPM(AliRsnMiniAnalysisTask *task, 
		Bool_t                 isMC, 
		AliPIDResponse::EBeamType collSys = AliPIDResponse::kPBPB, //=0, kPPB=1, kPBPB=2
		AliRsnCutSet           *cutsPair,             //cuts on the pair
		Bool_t                 enaMultSel = kTRUE,    //enable multiplicity axis
      		Float_t                masslow = 1.2,         //inv mass axis low edge 
		Float_t                massup = 3.2,          //inv mass axis upper edge 
		Int_t                  nbins = 2000,           //inv mass axis n bins
		Float_t                nsigma = 3.0,          //nsigma of TPC PID cut
		Bool_t                 enableMonitor = kTRUE) //enable single track QA plots
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
//  AliRsnCutPIDNSigma *iCutTPCTOFNSigma = new AliRsnCutPIDNSigma("CutTPCTOFNSigma", AliPID::kPion, AliRsnCutPIDNSigma::kTPC);
  //for setting PID cuts without any selection in pT
  //iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
 
//for setting PID cuts in given pT ranges
 iCutTPCNSigma->AddPIDRange(5.0, 0.0, 0.35);
 iCutTPCNSigma->AddPIDRange(3.0, 0.35, 0.5);
 iCutTPCNSigma->AddPIDRange(2.0, 0.5, 20.);

 iCutTOFNSigma->AddPIDRange(3.0, 0.0, 1.5);
 iCutTOFNSigma->AddPIDRange(2.5, 1.5, 1E20); 

//New, 23 april 2019
//        iCutTPCNSigma->SinglePIDRange(3.0);
//	//
//	iCutTOFNSigma->AddPIDRange(0.00,       0.00, 0.40);  
//	iCutTOFNSigma->AddPIDRange(3.0, 0.40, 1.e6);  
//	//
//        iCutTPCTOFNSigma->SinglePIDRange(5.0);

// end new
 //iCutTOFNSigma->SinglePIDRange(3.0);
 AliRsnCutSet * myCutSet = new AliRsnCutSet("MyCutSet", AliRsnTarget::kDaughter);
 myCutSet->AddCut(fCutQuality);
 myCutSet->AddCut(iCutTPCNSigma);
 myCutSet->AddCut(iCutTOFMatch);
 myCutSet->AddCut(iCutTOFNSigma);
// myCutSet->AddCut(iCutTPCTOFNSigma); 

// scheme: quality & [ (TPCsigma & !TOFmatch) | (TPCsigma & TOFsigma) ]
myCutSet->SetCutScheme( Form("%s & ((%s & (!%s))|(%s&%s))",fCutQuality->GetName(), iCutTPCNSigma->GetName(), iCutTOFMatch->GetName(),iCutTOFNSigma->GetName(), iCutTPCNSigma->GetName()) ) ;

//new cut scheme, 23 April 2019
// quality & [ ( TOFmatch & TOF & TPCTOF ) || ( TPConly ) ]
//  myCutSet->SetCutScheme( Form(" %s & ( ( %s & %s & %s ) | ( %s ) )",
//			 fCutQuality->GetName(),
//			 iCutTOFMatch->GetName(), iCutTOFNSigma->GetName(), iCutTPCTOFNSigma->GetName(),
//iCutTPCNSigma->GetName()) ) ;
// end!
 
//equivalent to
//SetCutScheme("CutQuality & ( (CutTPCNSigma & ! CutTOFMatch) | (CutTPCNSigma & CutTOFNSigma) ) ");
  Int_t icutPi = task->AddTrackCuts(myCutSet);
 
  //set daughter cuts
  Int_t icut1 = icutPi;
  //Int_t icut2 = icutPi;

  //monitor single-track selection based on track quality cuts only
//  AliRsnCutSetDaughterParticle * cutSetQuality = new AliRsnCutSetDaughterParticle("cutQuality", AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, 
//10.0, aodFilterBit, useCrossedRows);
//  Int_t icutQuality = task->AddTrackCuts(cutSetQuality);

   // selections for daugthers of Lambdas

  Float_t DCAxy;
  Float_t crossedRows = 70;
  Float_t rowsbycluster = 0.8;
  Float_t DCAxy = 0.15; //0.06;
     
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("qualityDaughterLs");
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.15,1.E10);
  esdTrackCuts->SetRequireTPCRefit();
  esdTrackCuts->SetAcceptKinkDaughters(0); //
  //esdTrackCuts->SetMinNCrossedRowsTPC(useCrossedRows);
  //esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(rowsbycluster);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMinDCAToVertexXY(DCAxy); //Use one of the two - pt dependent or fixed value cut.
  
  //V0s
  Float_t pi_Ls_PIDCut = 5.0;
  Float_t     LsDCA = 0.3;
  Float_t     LsCosPoinAn = 0.97;
  Float_t     LsDaughDCA=1.6;
  Float_t     massTol = 0.01;
  Float_t massTolVeto = 0.0043;
  Bool_t Switch = kFALSE;
  Float_t pLife = 50;
  Float_t v0rapidity= 0.5;
  Float_t  radiuslow=1.4;
  Int_t       tol_switch = 1;
  Double_t tol_sigma = 6;

  // selections for Lambdas
    AliRsnCutV0 *cutLambdas = new AliRsnCutV0("cutLambdas", kLambda0, AliPID::kProton, AliPID::kPion);
    cutLambdas->SetPIDCutPion(pi_Ls_PIDCut);        // PID for the pion daughter of Ls
    cutLambdas->SetPIDCutProton(pi_Ls_PIDCut);        // PID for the proton daughter of Ls
    cutLambdas->SetESDtrackCuts(esdTrackCuts);
    cutLambdas->SetMaxDaughtersDCA(LsDaughDCA);
    cutLambdas->SetMaxDCAVertex(LsDCA);
    cutLambdas->SetMinCosPointingAngle(LsCosPoinAn);
    cutLambdas->SetTolerance(massTol);
    //cutLambdas->SetToleranceVeto(massTolVeto);   //Rejection range for Competing V0 Rejection
    cutLambdas->SetSwitch(Switch);
    cutLambdas->SetfLife(pLife);
    cutLambdas->SetfLowRadius(radiuslow);
    cutLambdas->SetfHighRadius(100);
    cutLambdas->SetMaxRapidity(v0rapidity);
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
    cutAntiLambdas->SetMaxDCAVertex(LsDCA);
    cutAntiLambdas->SetMinCosPointingAngle(LsCosPoinAn);
    cutAntiLambdas->SetTolerance(massTol);
    //cutLambdas->SetToleranceVeto(massTolVeto);   //Rejection range for Competing V0 Rejection
    cutAntiLambdas->SetSwitch(Switch);
    cutAntiLambdas->SetfLife(pLife);
    cutAntiLambdas->SetfLowRadius(radiuslow);
    cutAntiLambdas->SetfHighRadius(100);
    cutAntiLambdas->SetMaxRapidity(v0rapidity);
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
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality       */ Int_t multID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
  /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt, kFALSE);
  /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt, kFALSE);
  /* 1st daughter p   */ Int_t fdp    = task->CreateValue(AliRsnMiniValue::kFirstDaughterP, kFALSE);
  /* 2nd daughter p   */ Int_t sdp    = task->CreateValue(AliRsnMiniValue::kSecondDaughterP, kFALSE);
  
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
    out->AddAxis(ptID, 200, 0.0, 20.0); //default use mother pt
    //axis Z: multiplicity
      if (enaMultSel) out->AddAxis(multID, 100, 0.0, 100.0);
  }
 
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

  // if (isMC) {
  //   //TRUE RECO PAIRS - TEMPLATE FOR BG
  //   for (Int_t ibg = 0; ibg <7; ibg++) {
  //     AliRsnMiniOutput * outtempl = task->CreateOutput(Form("bg_%s", bgTemplate[ibg].Data()), output.Data(),"TRUE");
  //     outtempl->SetCutID(0, icut1);
  //     outtempl->SetCutID(1, icut2);
  //     outtempl->SetCharge(0, bgTemplateC1[0]);
  //     outtempl->SetCharge(1, bgTemplateC2[0]);
  //     outtempl->SetDaughter(0, bgID1[ibg]);
  //     outtempl->SetDaughter(1, bgID2[ibg]);
  //     outtempl->SetMotherPDG(bgTemplatePDG[ibg]);
  //     outtempl->SetMotherMass(bgTemplateM[ibg]);
  //     outtempl->SetPairCuts(cutsPair);
  //     // axis X: invmass 
  //     outtempl->AddAxis(imID, nbins, masslow, massup);
  //     //axis Y: mother pt
  //     outtempl->AddAxis(ptID, 200, 0.0, 20.0); //default use mother pt
  //     // axis Z: multrality-multiplicity
  //     if (enaMultSel) outtempl->AddAxis(multID, 100, 0.0, 100.0);
  //   }
  //   //TRUE RECO PAIRS - MASS
  //   AliRsnMiniOutput * outtrue = task->CreateOutput(Form("trueSigPM_%s", partname.Data()), output.Data(),"TRUE");
  //   outtrue->SetCutID(0, icut1);
  //   outtrue->SetCutID(1, icut2);
  //   outtrue->SetCharge(0, charge1[0]);
  //   outtrue->SetCharge(1, charge2[0]);
  //   outtrue->SetDaughter(0, d1);
  //   outtrue->SetDaughter(1, d2);
  //   outtrue->SetMotherPDG(pdgCode);
  //   outtrue->SetMotherMass(mass);
  //   outtrue->SetPairCuts(cutsPair);
  //   // axis X: invmass 
  //   outtrue->AddAxis(imID, nbins, masslow, massup);
  //   //axis Y: mother pt
  //   outtrue->AddAxis(ptID, 200, 0.0, 20.0); //default use mother pt
  //   // axis Z: multiplicity
  //   if (enaMultSel) outtrue->AddAxis(multID, 100, 0.0, 100.0);

    
  //   //TRUE RECO PAIRS - MASS RESOLUTION
  //   AliRsnMiniOutput * outres = task->CreateOutput(Form("Mres_%s", partname.Data()), output.Data(),"TRUE");
  //   outres->SetCutID(0, icut1);
  //   outres->SetCutID(1, icut2);
  //   outres->SetCharge(0, charge1[0]);
  //   outres->SetCharge(1, charge2[0])
  //   outres->SetDaughter(0, d1);
  //   outres->SetDaughter(1, d2);
  //   outres->SetMotherPDG(pdgCode);
  //   outres->SetMotherMass(mass);
  //   outres->SetPairCuts(cutsPair);
  //   // axis X: invmass resolution
  //   outres->AddAxis(resID, 200, -0.01, 0.01);
  //   //axis Y: mother pt
  //   outres->AddAxis(ptID, 200, 0.0, 20.0);
  //   // axis Z: multiplicity
  //   if (enaMultSel) outres->AddAxis(multID, 100, 0.0, 100.0);

  //   //TRUE RECO PAIRS - rapidity
  //   AliRsnMiniOutput * outrap = task->CreateOutput(Form("trueRap_%s", partname.Data()), output.Data(),"TRUE");
  //   outrap->SetCutID(0, icut1);
  //   outrap->SetCutID(1, icut2);
  //   outrap->SetCharge(0, charge1[0]);
  //   outrap->SetCharge(1, charge2[0])
  //   outrap->SetDaughter(0, d1);
  //   outrap->SetDaughter(1, d2);
  //   outrap->SetMotherPDG(pdgCode);
  //   outrap->SetMotherMass(mass);
  //   outrap->SetPairCuts(cutsPair);
  //   outrap->AddAxis(ptID, 160, 0.0, 16.0);
  //   outrap->AddAxis(yID,  120, -0.6, 0.6);
  //   outrap->AddAxis(etaID, 200, -1., 1.);

    
  //   //GENERATED PAIRS
  //   AliRsnMiniOutput * outm = task->CreateOutput(Form("motherSigPM_%s", partname.Data()), output.Data(),"MOTHER");
  //   outm->SetDaughter(0, d1);
  //   outm->SetDaughter(1, d2);
  //   outm->SetMotherPDG(pdgCode);
  //   outm->SetMotherMass(mass);
  //   outm->SetPairCuts(cutsPair);
  //   outm->AddAxis(imID, nbins, masslow, massup);
  //   outm->AddAxis(ptID, 200, 0.0, 20.0);
  //   if (enaMultSel) outm->AddAxis(multID, 100, 0.0, 100.0);

  //   //GENERATED PAIRS
  //   AliRsnMiniOutput * outmy = task->CreateOutput(Form("motherRap_%s", partname.Data()), output.Data(),"MOTHER");
  //   outmy->SetDaughter(0, d1);
  //   outmy->SetDaughter(1, d2);
  //   outmy->SetMotherPDG(pdgCode);
  //   outmy->SetMotherMass(mass);
  //   outmy->SetPairCuts(cutsPair);
  //   outmy->AddAxis(ptID, 160, 0.0, 16.0);
  //   outmy->AddAxis(yID,  120, -0.6, 0.6);
  //   outmy->AddAxis(etaID, 200, -1., 1.);

  //   //f2 GENERATED PAIRS
  //   AliRsnMiniOutput * outm = task->CreateOutput("motherf2", output.Data(),"MOTHER");
  //   outm->SetDaughter(0, d1);
  //   outm->SetDaughter(1, d2);
  //   outm->SetMotherPDG(bgTemplatePDG[6]);
  //   outm->SetMotherMass(bgTemplateM[6]);
  //   outm->SetPairCuts(cutsPair);
  //   outm->AddAxis(imID, nbins, masslow, massup);
  //   outm->AddAxis(ptID, 200, 0.0, 20.0);
  //   if (enaMultSel) outm->AddAxis(multID, 100, 0.0, 100.0);

  //   //f2 GENERATED PAIRS
  //   AliRsnMiniOutput * outmy = task->CreateOutput("motherf2Rap", output.Data(),"MOTHER");
  //   outmy->SetDaughter(0, d1);
  //   outmy->SetDaughter(1, d2);
  //   outmy->SetMotherPDG(bgTemplatePDG[6]);
  //   outmy->SetMotherMass(bgTemplateM[6]);
  //   outmy->SetPairCuts(cutsPair);
  //   outmy->AddAxis(ptID, 160, 0.0, 16.0);
  //   outmy->AddAxis(yID,  120, -0.6, 0.6);
  //   outmy->AddAxis(etaID, 200, -1., 1.);
  // }

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
  // axisMass->SetBins(1.08,1.16,0.001);
axisMass->SetBins(1.08,3.16,0.001);
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

