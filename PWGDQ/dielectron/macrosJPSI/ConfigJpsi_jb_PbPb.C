void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);
void InitHF(AliDielectron* die, Int_t cutDefinition);

void SetupEventCuts(AliDielectron *die, ULong64_t triggers, Int_t cutDefinition);
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupV0Cuts( AliDielectron *die,  Int_t cutDefinition);
void SetupPairCuts( AliDielectron *die,  Int_t cutDefinition);

void ConfigEvtPlane(AliDielectron *die,  Int_t cutDefinition);
void ConfigBgrd(    AliDielectron *die,  Int_t cutDefinition);

void AddMCSignals(AliDielectron *die);
void SetEtaCorrection();
TVectorD *GetRunNumbers();
TVectorD *GetDeltaPhiBins();

TString names=("TPC;TOF;TRD;TOFTRD;TOFTRD2D;rec;Gam0;Gam01;Gam05;Gam10;Gam15;Gam20;EtaGap01;EtaGap02;EtaGap03;EtaGap04;EtaGap05;SubLS;SubRndm");
enum { kTPC=0, kTOF, kTRD, kTOFTRD, kTOFTRD2D, krec,kGam0, kGam01, kGam05, kGam10, kGam15, kGam20, kEtaGap01, kEtaGap02, kEtaGap03, kEtaGap04, kEtaGap05, kSubLS, kSubRndm };

TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntries();

Bool_t  isESD = kTRUE;
Bool_t  hasMC = kFALSE;
TString list  = gSystem->Getenv("LIST");

AliDielectron* ConfigJpsi_jb_PbPb(Int_t cutDefinition, TString prod="", ULong64_t triggers=AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB)
{
  //
  // Setup the instance of AliDielectron
  //

  // gsi train?
  TString trainRoot = gSystem->Getenv("TRAIN_ROOT");
  Bool_t isGSItrain = (trainRoot.IsNull()?kFALSE:kTRUE); 

  // find mc or not?
  if( list.IsNull()) list=prod;
  if( list.Contains("LHC10h")   || list.Contains("LHC11h")   ) hasMC=kFALSE;
  if( list.Contains("LHC11a10") || list.Contains("LHC12a17") ) hasMC=kTRUE;

  //ESD handler?
  isESD=(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliESDInputHandler::Class());

  // switch configurations ON and OFF
  if(hasMC) { // MONTE CARLO
    switch(cutDefinition) {
    case kTOF: /* */ break;
    default:         return 0x0;
    }
  } else { // COLLISION DATA
    switch(cutDefinition) {
    case kTOF: /* */ break;
    default:         return 0x0;
    }
  }

  // task name
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast())  name=arrNames->At(cutDefinition)->GetName();
  printf(" Adding %s%s config %s for %s \n",(isESD?"ESD":"AOD"),(hasMC?" MC":""),name.Data(),list.Data());

  // init AliDielectron
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("Track cuts: %s",name.Data()));
  die->SetHasMC(hasMC);

  // cut setup
  SetupEventCuts(die,triggers,cutDefinition);
  SetupTrackCuts(die,cutDefinition);
  SetupV0Cuts(die,cutDefinition);
  //  SetupPairCuts(die,cutDefinition);

  // Monte Carlo Signals
  if(hasMC) {
    AddMCSignals(die);
    printf(" Add %d MC signals \n",die->GetMCSignals()->GetEntriesFast());
  }

  // TRD efficiencies tables
  /*
    if (hasMC && list.Contains("LHC11a") ) {
     TString pidTab="$TRAIN_ROOT/util/dielectron/dielectron/TRDpidEff_eleProb07_TRDntr4_6.root";
     if(!isGSItrain) pidTab="$ALICE_ROOT/PWGDQ/dielectron/files/TRDpidEff_eleProb07_TRDntr4_6.root";

     if (gSystem->AccessPathName(gSystem->ExpandPathName(pidTab.Data())))
     Error("ConfigPbPb","PID table not found: %s",pidTab.Data());
     else
     die->SetTRDcorrectionFilename(pidTab.Data());
     }
  */

  // histogram setup
  InitHistograms(die,cutDefinition);
  printf(" Add %d class types to the histo manager \n",die->GetHistogramList()->GetEntries());

  // HF array setup
  InitHF(die,cutDefinition);

  // CF container setup, switched off
  InitCF(die,cutDefinition);
  /* printf(" Add %d pair, %d leg vars, %p steps and %p bins to the container \n",
		 die->GetCFManagerPair()->GetNvarsPair(),die->GetCFManagerPair()->GetNvarsLeg(),
		 die->GetCFManagerPair()->GetContainer(), die->GetCFManagerPair()->GetContainer() );
  */

  // bgrd estimators
  if(!hasMC) ConfigBgrd(die,cutDefinition);

  // tpc event plane configuration
  ConfigEvtPlane(die,cutDefinition);

  // prefilter settings
  die->SetPreFilterUnlikeOnly();
  //die->SetPreFilterAllSigns();
  //die->SetNoPairing();

  // setup eta correction
  // if(isESD && list.Contains("LHC10h")) SetEtaCorrection();

  // VZERO calibration
  if (isGSItrain && list.Contains("LHC10h")) {
    die->SetVZEROCalibrationFilename("$TRAIN_ROOT/util/dielectron/dielectron/VzeroCalibrationLHC10h.root");
    die->SetVZERORecenteringFilename("$TRAIN_ROOT/util/dielectron/dielectron/VzeroRecenteringLHC10h.root");
  }

  return die;
}

//______________________________________________________________________________________
void SetupEventCuts(AliDielectron *die, ULong64_t triggers, Int_t cutDefinition)
{
  //
  // Setup the event cuts
  //

  // trigger specific centrality cuts (reject trigger inefficiencies)
  Double_t minCent=0.0, maxCent=100.;
  if(!hasMC) {
    switch(triggers) {
    case AliVEvent::kCentral:     minCent= 0.; maxCent= 9.; break;
    case AliVEvent::kSemiCentral: minCent=12.; maxCent=53.; break;
    case AliVEvent::kMB:          minCent= 0.; maxCent=80.; break;
    default:                      minCent= 0.; maxCent=80.; break;
    }
  }
  //  if(cutDefinition >= kEtaGap01) {minCent=20.; maxCent=50.;} // v2 analysis

  // VZERO multiplicity vs. number ob global tracks cut
  TF1 *fMean  = new TF1("fMean", "pol1",               0,25e+3);
  fMean->SetParameters(691.633, 1.4892);
  TF1 *fSigma = new TF1("fSigma","[0]+sqrt([1]*x+[2])",0,25e+3);
  fSigma->SetParameters(-83.6599, 36.7677, 69530.7);

  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","eventCuts");
  if(!isESD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,+10.);
  eventCuts->SetCentralityRange(minCent,maxCent);
  eventCuts->SetCutOnV0MultipicityNTrks(fMean, fSigma, 4.0);
  eventCuts->Print();
  die->GetEventFilter().AddCuts(eventCuts);

}

//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //

  // Quality cuts
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  die->GetTrackFilter().AddCuts(cuts);

  // AOD track filter (needs to be first cut to speed up)
  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");
  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  //  trkFilter->SetMinNCrossedRowsOverFindable(0.6);
  if(!isESD) cuts->AddCut(trkFilter);

  //Pt cut, should make execution a bit faster
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("PtCut","PtCut");
  pt->AddCut(AliDielectronVarManager::kPt,1.1,1e30);    //0.8
  cuts->AddCut(pt);

  // track cuts ESD and AOD
  AliDielectronVarCuts *varCuts = new AliDielectronVarCuts("VarCuts","VarCuts");
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.9,   0.9); // -0.9, 0.9
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,     70.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.0);
  //  varCuts->AddCut(AliDielectronVarManager::kV0Index0,     0.0);
  if(cutDefinition==kTOF || cutDefinition==kTOFTRD || cutDefinition==kTOFTRD2D)
    varCuts->AddCut(AliDielectronVarManager::kTOFbeta,      0.2,   0.9, kTRUE);
  cuts->AddCut(varCuts);
  varCuts->Print();

  AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  trkCuts->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 15); // ITS-4 = 1+2+4+8
  // trkCuts->SetMaxWaivedITSNcls(1);
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE);
  cuts->AddCut(trkCuts);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronPID *pid = new AliDielectronPID("PID","PID");
  ////////////////////////////////// DATA
  if(!hasMC) {
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,4.0,0.,0.,kTRUE);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-100.,3.5,0.,0.,kTRUE);

    if(cutDefinition==kTRD || cutDefinition==kTOFTRD)
      pid->AddCut(AliDielectronPID::kTRDeleEff,AliPID::kElectron,.8,1.,3.5.,6.,kFALSE,
                  AliDielectronPID::kIfAvailable,AliDielectronVarManager::kTRDpidQuality);

    if(cutDefinition==kTOFTRD2D)
      pid->AddCut(AliDielectronPID::kTRDeleEff2D,AliPID::kElectron,.8,1.,3.5.,6.,kFALSE,
                  AliDielectronPID::kIfAvailable,AliDielectronVarManager::kTRDpidQuality);
  }

  ////////////////////////////////// MC
  if(hasMC) {

    // electron
    Double_t nSigmaPi = 3.5; Double_t nSigmaP = 3.5;
    Double_t resolution=0.0549;
    Double_t BBpro[5] = {0};
    Double_t BBpio[5] = {0};

    for(Int_t icent=0; icent<8; icent++) {

      switch (icent) {
        case 0:  // 0-10%
          BBpro[0] = 0.031555;  BBpro[1] = 26.0595; BBpro[2] = 3.02422e-11;  BBpro[3] = 2.05594; BBpro[4] = 5.99848;
          BBpio[0] = 0.0252122; BBpio[1] = 38.8991; BBpio[2] = 4.0901e-11;   BBpio[3] = 5.27988; BBpio[4] = 4.3108;
          break;
        case 1:  // 10-20%
          BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
          BBpio[0] = 0.0252127; BBpio[1] = 33.8617; BBpio[2] = 3.56866e-11;  BBpio[3] = 5.24831; BBpio[4] = 4.31093;
          break;
        case 2:  // 20-30%
          BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
          BBpio[0] = 0.0263205; BBpio[1] = 37.9307; BBpio[2] = 4.29724e-11;  BBpio[3] = 5.74458; BBpio[4] = 4.32459;
          break;
        case 3:  // 30-40%
          BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
          BBpio[0] = 0.026294;  BBpio[1] = 39.0346; BBpio[2] = 4.12261e-11;  BBpio[3] = 5.28808; BBpio[4] = 4.31301;
          break;
        case 4:  // 40-50%
          BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
          BBpio[0] = 0.0263134; BBpio[1] = 38.2084; BBpio[2] = 3.75159e-11;  BBpio[3] = 5.78125; BBpio[4] = 4.31363;
          break;
        case 5:  // 50-60%
          BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
          BBpio[0] = 0.0263134; BBpio[1] = 38.2084; BBpio[2] = 3.75159e-11;  BBpio[3] = 5.78125; BBpio[4] = 4.31363;
          break;
        case 6:  // 60-70%
          BBpro[0] = 0.031555;  BBpro[1] = 26.0595; BBpro[2] = 3.02422e-11;  BBpro[3] = 2.05594; BBpro[4] = 5.99848;
          BBpio[0] = 0.026302;  BBpio[1] = 38.6888; BBpio[2] = 3.56792e-11;  BBpio[3] = 5.2465;  BBpio[4] = 4.31094;
          break;
        case 7:  // 70-80%
          BBpro[0] = 0.0315171; BBpro[1] = 25.8656; BBpro[2] = 3.03896e-11;  BBpro[3] = 2.05802; BBpro[4] = 5.99999;
          BBpio[0] = 0.0263134; BBpio[1] = 38.2084; BBpio[2] = 3.75159e-11;  BBpio[3] = 5.78125; BBpio[4] = 4.31363;
          break;
        case 8:  // 80-90%
          BBpro[0] = 0.0313438; BBpro[1] = 25.8666; BBpro[2] = 4.5457e-11;   BBpro[3] = 2.07912; BBpro[4] = 5.99986;
          BBpio[0] = 0.0252127; BBpio[1] = 33.8617; BBpio[2] = 3.56866e-11;  BBpio[3] = 5.24831; BBpio[4] = 4.31093;
          break;
        case 9:  // 90-100%
          BBpro[0] = 0.0319126; BBpro[1] = 36.8784; BBpro[2] = 3.4274e-11;   BBpro[3] = 3.2431;  BBpro[4] = 5.93388;
          BBpio[0] = 0.027079;  BBpio[1] = 67.5936; BBpio[2] = 9.72548e-11;  BBpio[3] = 9.61382; BBpio[4] = 5.99372;
          break;
      }


      TF1 *ffPro=new TF1(Form("fBethe%d_c%d",AliPID::kProton,icent), Form("(%f*%f+(AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])-AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])))/%f", nSigmaP,resolution, AliPID::ParticleMass(AliPID::kProton), AliPID::ParticleMass(AliPID::kElectron), resolution), 0.05,200.);

      TF1 *ffPio=new TF1(Form("fBethe%d_c%d",AliPID::kPion,icent), Form("(%f*%f+(AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])-AliExternalTrackParam::BetheBlochAleph(x/%f,[0],[1],[2],[3],[4])))/%f", nSigmaPi,resolution, AliPID::ParticleMass(AliPID::kPion), AliPID::ParticleMass(AliPID::kElectron), resolution), 0.05,200.);

      //LHC11a10b
      if (list.Contains("LHC11a10b")) {
        printf(" LHC11a10b parameters\n");
        ffPro->SetParameters(BBpro[0],BBpro[1],BBpro[2],BBpro[3],BBpro[4]);
        ffPio->SetParameters(BBpio[0],BBpio[1],BBpio[2],BBpio[3],BBpio[4]);

        // proton cut
        pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,ffPro,10,((double)icent)*10.,((double)icent+1)*10,
                    kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kCentrality);
        // pion cut
        pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,ffPio,10,((double)icent)*10.,((double)icent+1)*10,
                    kFALSE,AliDielectronPID::kRequire,AliDielectronVarManager::kCentrality);
      }
    }

    // shifts for the nSigma electrons
    TGraph* nSigmaCorrection = new TGraph();
    // LHC11a10b
    if (list.Contains("LHC11a10b")) {
      nSigmaCorrection->SetPoint(0, 137161., -0.50-(0.28));
      nSigmaCorrection->SetPoint(1, 139510., -0.50-(0.28));
      pid->SetCorrGraph(nSigmaCorrection);
    }

  } //hasMC

  ////////////////////////////////// DATA + MC
  // pid cuts TPC + TOF & TRD
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-4.,4.);

  if(cutDefinition==kTOF || cutDefinition==kTOFTRD || cutDefinition==kTOFTRD2D)
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3,3.,0.,0.,kFALSE,AliDielectronPID::kIfAvailable);

  if(cutDefinition!=krec) cuts->AddCut(pid);
  /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ PID CUTS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/


  // exclude conversion electrons selected by the tender
  AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
  noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
  //  cuts->AddCut(noconv);

}

//______________________________________________________________________________________
void SetupV0Cuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the V0 cuts
  //

  // Quality cuts
  TIter next(die->GetTrackFilter().GetCuts());
  AliAnalysisCuts *cuts;
  while((cuts = (AliAnalysisCuts*)next())) {
    if(cuts->IsA() == AliDielectronCutGroup::Class())  break;
  }

  AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("IsGamma","IsGamma");
  gammaV0Cuts->SetPdgCodes(22,11,11);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE);
  //  gammaV0Cuts->AddCut(AliDielectronVarManager::kOpeningAngle,              0.0,   0.1, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35,  0.35, kFALSE);
  gammaV0Cuts->SetExcludeTracks(kTRUE);
  gammaV0Cuts->Print();

 //  const Double_t |cutAlphaG| < 0.35; &&  const Double_t cutQTG < 0.05;                                                        
 //  const Double_t |cutAlphaG2|[2] = {0.6, 0.8}; &&  const Double_t cutQTG2 < 0.04;  

  if(cuts)
    ((AliDielectronCutGroup*)cuts)->AddCut(gammaV0Cuts);
  else
    die->GetTrackFilter().AddCuts(gammaV0Cuts);
}

//______________________________________________________________________________________
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the pair cuts
  //

  // conversion rejection
  Double_t gCut;
  switch(cutDefinition) {
    case kTPC:    gCut=0.05;  break;
    case krec:    gCut=0.05;  break;
    case kGam10:  gCut=0.1;   break;
    case kGam15:  gCut=0.15;  break;
    case kGam20:  gCut=0.2;   break;
    case kGam05:  gCut=0.05;  break;
    case kGam01:  gCut=0.01;  break;
    case kGam0:   gCut=0.0;   break;
    default: gCut=0.05;             // default
  }

  AliDielectronVarCuts *gammaCuts = new AliDielectronVarCuts("GammaCuts","GammaCuts");
//  gammaCuts->AddCut(AliDielectronVarManager::kOpeningAngle, 0.0,   0.1,  kTRUE);
//  gammaCuts->AddCut(AliDielectronVarManager::kLegDist,      0.0,   0.25, kTRUE);
//  gammaCuts->AddCut(AliDielectronVarManager::kR,            3.0,   90.0, kTRUE);
//  gammaCuts->AddCut(AliDielectronVarManager::kPsiPair,      0.0,   0.05, kTRUE);
//  gammaCuts->AddCut(AliDielectronVarManager::kChi2NDF,      0.0,   10.0, kTRUE);
  gammaCuts->AddCut(AliDielectronVarManager::kM,            0.0,   gCut);
  die->GetPairPreFilter().AddCuts(gammaCuts);


  // rapidity selection
  //  AliDielectronVarCuts *rapCut=new AliDielectronVarCuts("|Y|<.9","|Y|<.9");
  // rapCut->AddCut(AliDielectronVarManager::kY,-0.9,0.9);
  // die->GetPairFilter().AddCuts(rapCut);

}

//______________________________________________________________________________________
void ConfigBgrd(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Configurate the background estimators
  //

  // add track rotations
  AliDielectronTrackRotator *rot=new AliDielectronTrackRotator;
  rot->SetIterations(10);
  rot->SetConeAnglePhi(TMath::Pi());
  rot->SetStartAnglePhi(TMath::Pi());
  //      die->SetTrackRotator(rot);

  // add mixed events
  AliDielectronMixingHandler *mix=new AliDielectronMixingHandler;
  mix->AddVariable(AliDielectronVarManager::kZvPrim,      "-10.,-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.,10.");
  //  mix->AddVariable(AliDielectronVarManager::kCentrality,  8,  0.,80.);
  mix->AddVariable(AliDielectronVarManager::kTPCrpH2,     8,  TMath::Pi()/-2., TMath::Pi()/2.);
  mix->AddVariable(AliDielectronVarManager::kTPCmagH2,    "0.,20.,50.,80.,110.,150.,500.");
  mix->SetMixType(AliDielectronMixingHandler::kAll);
  mix->SetDepth(150);
  die->SetMixingHandler(mix);

}

//______________________________________________________________________________________
void ConfigEvtPlane(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Configurate the TPC event plane 
  //

  Double_t gGap;
  switch(cutDefinition) {
  case kEtaGap01:   gGap=0.1;   break;
  case kEtaGap02:   gGap=0.2;   break;
  case kEtaGap03:   gGap=0.3;   break;
  case kEtaGap04:   gGap=0.4;   break;
  case kEtaGap05:   gGap=0.5;   break;
  default: gGap=0.0;
  }

  AliDielectronVarCuts *poi = new AliDielectronVarCuts("PoI","PoI");
  poi->AddCut(AliDielectronVarManager::kM,2.92,3.20);     // particles of interest, jpsi mass window
  die->GetEventPlanePOIPreFilter().AddCuts(poi);

  // eta gap in tpc event plane
  if(cutDefinition >= kEtaGap01 &&
     cutDefinition <  kSubLS     ) {
    AliDielectronVarCuts *etaGap = new AliDielectronVarCuts(AliDielectronVarManager::GetValueName(AliDielectronVarManager::kEta),"etaGap");
    etaGap->AddCut(AliDielectronVarManager::kEta,-1*gGap,gGap,kTRUE);
    die->GetEventPlanePreFilter().AddCuts(etaGap);
  }

  if(cutDefinition==kSubLS) die->SetLikeSignSubEvents();
  die->SetPreFilterEventPlane();
}

//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //

  // booleans for histo selection
  Bool_t bHistTrackQA=kFALSE, bHistEvts = kFALSE, bHistPair = kFALSE, bHistPairME = kFALSE, bHistFlow = kFALSE, bHistFlowQA=kFALSE;
  switch (cutDefinition) {
  case kTPC:
  case kTOF:      bHistEvts=kTRUE; bHistFlow=kTRUE; bHistPair=kTRUE; bHistPairME=kFALSE; break;
  case kTRD:
  case kTOFTRD:
  case kTOFTRD2D: bHistFlow=kTRUE; bHistPair=kTRUE; break;
  case krec:      break;
    //  case kITScls:
    //  case kITSamy:
    //  case kDCA:
    //  case kChi:      break;
  case kEtaGap01:
  case kEtaGap02:
  case kEtaGap03:
  case kEtaGap04:
  case kEtaGap05:
  case kSubRndm:  bHistFlow=kTRUE; bHistFlowQA=kTRUE; break;
  case kSubLS:    bHistFlow=kTRUE; bHistFlowQA=kTRUE; bHistEvts=kTRUE; break;
  case kGam0:
  case kGam01:
  case kGam05:
  case kGam10:
  case kGam15:
  case kGam20:    break;
  }


  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(die->GetName(),die->GetTitle());

  //add histograms to event class
  histos->AddClass("Event");
  histos->UserHistogram("Event","","", 100, 0.0, 100.0,   AliDielectronVarManager::kCentrality);

  ////// EVENT HISTOS /////
  if(bHistEvts) {
    histos->UserHistogram("Event","","", GetRunNumbers(), AliDielectronVarManager::kRunNumber);
    histos->UserHistogram("Event","","", 300,-15.,15.,	  AliDielectronVarManager::kZvPrim);
    histos->UserProfile(  "Event","","", AliDielectronVarManager::kNacc,        80, 0., 80.,	 AliDielectronVarManager::kCentrality);
    histos->UserProfile(  "Event","","", AliDielectronVarManager::kNVtxContrib, 80, 0., 80.,   AliDielectronVarManager::kCentrality);
  } //hist: event


  ////// FLOW //////
  if(bHistFlow) {
    // EP Qvector magnitudes // TODO move to QA
    histos->UserHistogram("Event","","", 200,0.,200., AliDielectronVarManager::kTPCmagH2uc);
    histos->UserHistogram("Event","","", 200,0.,800., AliDielectronVarManager::kv0ACmagH2);
    histos->UserHistogram("Event","","", 200,0.,800., AliDielectronVarManager::kv0AmagH2);
    histos->UserHistogram("Event","","", 200,0.,800., AliDielectronVarManager::kv0CmagH2);
    // RP angles versus centrality
    histos->UserHistogram("Event","","", 16,0.,80.,100,-2.,2., AliDielectronVarManager::kCentrality,AliDielectronVarManager::kTPCrpH2);
    histos->UserHistogram("Event","","", 16,0.,80.,100,-2.,2., AliDielectronVarManager::kCentrality,AliDielectronVarManager::kTPCsub1rpH2);
    histos->UserHistogram("Event","","", 16,0.,80.,100,-2.,2., AliDielectronVarManager::kCentrality,AliDielectronVarManager::kTPCsub2rpH2);
    histos->UserHistogram("Event","","", 16,0.,80.,100,-2.,2., AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0ArpH2);
    histos->UserHistogram("Event","","", 16,0.,80.,100,-2.,2., AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0CrpH2);
    histos->UserHistogram("Event","","", 16,0.,80.,100,-2.,2., AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0ACrpH2);
    histos->UserHistogram("Event","","", 16,0.,80.,100,-2.,2., AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0A0rpH2);
    histos->UserHistogram("Event","","", 16,0.,80.,100,-2.,2., AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0C0rpH2);
    histos->UserHistogram("Event","","", 16,0.,80.,100,-2.,2., AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0A3rpH2);
    histos->UserHistogram("Event","","", 16,0.,80.,100,-2.,2., AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0C3rpH2);
  } // hist: flow

  if(bHistFlowQA) {
    // TPC event plane
    histos->UserHistogram("Event","","", 100,-1500.,1500., AliDielectronVarManager::kTPCxH2);
    histos->UserHistogram("Event","","", 100,-1500.,1500., AliDielectronVarManager::kTPCyH2);
    histos->UserHistogram("Event","","", 100,   -2.,   2., AliDielectronVarManager::kTPCrpH2);
    histos->UserHistogram("Event","","", 100,   -2.,   2., AliDielectronVarManager::kTPCsub1rpH2);
    histos->UserHistogram("Event","","", 100,   -2.,   2., AliDielectronVarManager::kTPCsub2rpH2);
    histos->UserHistogram("Event","","", 100,   -1.,   1., AliDielectronVarManager::kTPCsub12DiffH2);
    // EP resolution calculation
    histos->UserHistogram("Event","","", 80,0.,80., 300,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0ATPCDiffH2);
    histos->UserHistogram("Event","","", 80,0.,80., 300,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0CTPCDiffH2);
    histos->UserHistogram("Event","","", 80,0.,80., 300,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0Av0CDiffH2);
    histos->UserHistogram("Event","","", 80,0.,80., 300,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kTPCsub12DiffH2);
    histos->UserHistogram("Event","","", 80,0.,80., 300,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0Av0C0DiffH2);
    histos->UserHistogram("Event","","", 80,0.,80., 300,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0Av0C3DiffH2);
    histos->UserHistogram("Event","","", 80,0.,80., 300,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0Cv0A0DiffH2);
    histos->UserHistogram("Event","","", 80,0.,80., 300,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0Cv0A3DiffH2);
    histos->UserHistogram("Event","","", 80,0.,80., 300,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0A0v0A3DiffH2);
    histos->UserHistogram("Event","","", 80,0.,80., 300,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kv0C0v0C3DiffH2);
    // detector effects
    histos->UserHistogram("Event","","", 10,0.,100., 300,-1.0,1.0, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kTPCsub12DiffH2Sin);
    // recentering stuff
    histos->UserProfile("Event","","", AliDielectronVarManager::kTPCxH2,
			AliDielectronHelper::MakeLinBinning(8, 0.,80.), GetRunNumbers(),
			AliDielectronVarManager::kCentrality, AliDielectronVarManager::kRunNumber);
    histos->UserProfile("Event","","", AliDielectronVarManager::kTPCyH2,
			AliDielectronHelper::MakeLinBinning(8, 0.,80.), GetRunNumbers(),
			AliDielectronVarManager::kCentrality, AliDielectronVarManager::kRunNumber);
  } //hist: flowQA

  if(bHistPair) {
    //Initialise histogram classes
    histos->SetReservedWords("Track;Pair");

    //Pair classes
    // to fill also mixed event histograms loop until 7 or 10
    for (Int_t i=0; i<(bHistPairME ? 8 : 3); ++i){
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
    }
    //legs from pair (fill SE)
    for (Int_t i=0; i<3; ++i){
      histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i)));
    }
    //Track classes
    //to fill also track info from 2nd event loop until 2
    for (Int_t i=0; i<2; ++i){
      histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
    }
    //track rotation
    //   histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
    //   histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));

    // Vertex
    histos->UserHistogram("Track","","", 500,-1.,1., AliDielectronVarManager::kImpactParXY);
    histos->UserHistogram("Track","","", 600,-3.,3., AliDielectronVarManager::kImpactParZ);
    // Kinematics
    histos->UserHistogram("Track","","", 400,0,20.,  AliDielectronVarManager::kPt);
    histos->UserHistogram("Track","","", 200,-1,1, 200,0,6.285, AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
    // TPC
    histos->UserHistogram("Track","","", 160,-0.5,159.5, AliDielectronVarManager::kNclsTPC);
    histos->UserHistogram("Track","","", 160,-0.5,159.5, AliDielectronVarManager::kTPCsignalN);
    histos->UserHistogram("Track","","", 160,-0.5,159.5, AliDielectronVarManager::kNFclsTPCr);
    histos->UserHistogram("Track","","", 160,-0.5,159.5, 160,-0.5,159.5, AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
    // TRD
    histos->UserHistogram("Track","","",   8,-0.5,  7.5, AliDielectronVarManager::kTRDpidQuality);
    // PID
    histos->UserHistogram("Track","","", 400,0.2,20.,200,0.,200., AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);
    histos->UserHistogram("Track","","", 400,0.2,20.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
    histos->UserHistogram("Track","","", 250,0.0,5.0,300,0.,1.2,  AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta,kTRUE);

    histos->UserHistogram("Pair","","", 210,-1.05,1.05, 100,0.,2.5, AliDielectronVarManager::kArmAlpha,AliDielectronVarManager::kArmPt);

    ///// add histograms to Pair classes /////
    histos->UserHistogram("Pair","","",  300,.0,300*0.04, AliDielectronVarManager::kM); // 40MeV bins, 12GeV/c2
    histos->UserHistogram("Pair","","",  100,-1.,1.,      AliDielectronVarManager::kY);
    histos->UserHistogram("Pair","","",  400,0,20.,       AliDielectronVarManager::kPt);
    histos->UserHistogram("Pair","","",  100,0.,3.15,     AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","","",  100,0.,20,       AliDielectronVarManager::kChi2NDF);
    histos->UserHistogram("Pair","","",  100,0.,3.15,     AliDielectronVarManager::kPsiPair);
    histos->UserHistogram("Pair","","",  200,0.,100.,     AliDielectronVarManager::kR);
    histos->UserHistogram("Pair","","",   50,0.,5.,       AliDielectronVarManager::kLegDist);
    histos->UserHistogram("Pair","","",   50,0.,5.,       AliDielectronVarManager::kLegDistXY);

    //// FLOW results use tprofiles
    if(bHistFlow) {
      histos->UserProfile("Pair","",
			  "cos(2(#varphi-#Psi^{V0AC}));mass (GeV/c^{2});centrality (%);p_{T} (GeV/c)",
			  AliDielectronVarManager::kv0ACrpH2FlowV2,
			  125,0.,125*.04, 10, 0.,100., 200,0.,100.,
			  AliDielectronVarManager::kM, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kPt);
      histos->UserProfile("Pair","",
			  "cos(2(#varphi-#Psi^{V0A}));mass (GeV/c^{2});centrality (%);p_{T} (GeV/c)",
			  AliDielectronVarManager::kv0ArpH2FlowV2,
			  125,0.,125*.04, 10, 0.,100., 200,0.,100.,
			  AliDielectronVarManager::kM, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kPt);
      histos->UserProfile("Pair","",
			  "cos(2(#varphi-#Psi^{V0C}));mass (GeV/c^{2});centrality (%);p_{T} (GeV/c)",
			  AliDielectronVarManager::kv0CrpH2FlowV2,
			  125,0.,125*.04, 10, 0.,100., 200,0.,100.,
			  AliDielectronVarManager::kM, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kPt);

      histos->UserProfile("Pair","","",			  AliDielectronVarManager::kv0ArpH2FlowV2,
			  125,0.,125*.04, 10, 0.,100.,	  AliDielectronVarManager::kM, AliDielectronVarManager::kCentrality);

    } //hist: flow
  } //hist: pair


  ////// MONTE CARLO //////
  /*
    if(cutDefinition == kTOFTRD && hasMC) {
    histos->AddClass("MCEvent");
    histos->UserHistogram("MCEvent","Cent_NJPsis","Centrality vs. generated incl. J/#psi per event;centrality (%);N_{J/#psi}",
    10,0.,100., 21,-0.5,20.5,
    AliDielectronVarManager::kCentrality,AliDielectronVarManager::kNumberOfJPsis);
    }
  */

  //  histos->UserHistogram("Track","TRDprobEle",";P_{ele}^{TRD};#tracks", 
  //			100,0.,1.,            AliDielectronVarManager::kTRDprobEle);
  
  die->SetHistogramManager(histos);
}

void InitHF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setup the HF arrays
  //


  AliDielectronHF *hf=new AliDielectronHF(die->GetName(),die->GetTitle());
  //  if(hasMC) hf->SetStepForMCGenerated();
  hf->SetPairTypes(AliDielectronHF::kAll);
  hf->SetVariable(AliDielectronVarManager::kM, 125, 0.0, 0.04*125);

  hf->AddCutVariable(AliDielectronVarManager::kCentrality,  "0.,5.,10.,20.,40.,50.,60.,80."  );
  hf->AddCutVariable(AliDielectronVarManager::kPt,          "0.,2.5,5.,100."                 );
  hf->AddCutVariable(AliDielectronVarManager::kDeltaPhiv0ArpH2, 8,-1.*TMath::Pi(),TMath::Pi());
  hf->AddCutVariable(AliDielectronVarManager::kDeltaPhiv0CrpH2, 8,-1.*TMath::Pi(),TMath::Pi());
  //  hf->AddCutVariable(AliDielectronVarManager::kY,           1, -0.9, 0.9                     );
  //  hf->AddCutVariable(AliDielectronVarManager::kPt,          "0.8, 1.0, 1.1, 1.2, 1.5, 100.0", kTRUE, AliDielectronHF::kBinToMax);
  // hf->AddCutVariable(AliDielectronVarManager::kNclsTPC,     "70,90,100,120,160",              kTRUE, AliDielectronHF::kBinToMax);
  //  hf->AddCutVariable(AliDielectronVarManager::kTPCnSigmaEle,"-4,-3,-2.5,-2,2,2.5,3,4",        kTRUE, AliDielectronHF::kSymBin);
  //hf->AddCutVariable(AliDielectronVarManager::kTPCnSigmaPio,"3.,3.5,4.,100.",                 kTRUE, AliDielectronHF::kBinToMax);
  //hf->AddCutVariable(AliDielectronVarManager::kITSLayerFirstCls,4,0.,4.,              kFALSE, kTRUE, AliDielectronHF::kBinFromMin);
  //hf->AddCutVariable(AliDielectronVarManager::kNclsITS,         5,2.,7.,              kFALSE, kTRUE, AliDielectronHF::kBinToMax);
  //hf->AddCutVariable(AliDielectronVarManager::kRunNumber,  GetRunNumbers());

  die->SetHistogramArray(hf);
}

void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setup the CF Manager if needed
  //

  // switch off ep systematics
  if(cutDefinition >  kEtaGap01 && cutDefinition != kSubRndm  ) continue;

  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());

  // event variables
  cf->AddVariable(AliDielectronVarManager::kCentrality,      "0.,5.,10.,20.,40.,50.,60.,80.");
  // pair variables
  cf->AddVariable(AliDielectronVarManager::kPairType,        11,    0, 11);
  cf->AddVariable(AliDielectronVarManager::kM,               125, 0.0, 0.04*125); // 40MeV bins
  cf->AddVariable(AliDielectronVarManager::kPt,              "0., 1., 2.5, 5., 100.0");
  cf->AddVariable(AliDielectronVarManager::kDeltaPhiv0ArpH2, GetDeltaPhiBins());
  cf->AddVariable(AliDielectronVarManager::kDeltaPhiv0CrpH2, GetDeltaPhiBins());

  /*
  //    if(!hasMC) cf->AddVariable(AliDielectronVarManager::kZvPrim,20, -10., 10.);
  if(hasMC)  cf->AddVariable(AliDielectronVarManager::kNacc,20,0.,3000.0);
  if(hasMC)  cf->AddVariable(AliDielectronVarManager::kNVtxContrib,20,0.,4000.);
  if(hasMC)  cf->AddVariable(AliDielectronVarManager::kRunNumber, GetRunNumbers() );

  // pair variables
  //    cf->AddVariable(AliDielectronVarManager::kY,"-0.8,0.8");
  cf->AddVariable(AliDielectronVarManager::kM,125,0.,125*.04); //40Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,11,0,11);
  cf->AddVariable(AliDielectronVarManager::kPt,"0., 1., 2.5, 5., 100.0");
  if(hasMC) cf->AddVariable(AliDielectronVarManager::kTRDpidEffPair,101,0.0,1.01);
  //    if(hasMC) cf->AddVariable(AliDielectronVarManager::kThetaCS,15,-1.,1.);

  // flow variables
  cf->AddVariable(AliDielectronVarManager::kDeltaPhiv0ArpH2,4,-1.*TMath::Pi(),TMath::Pi());
  cf->AddVariable(AliDielectronVarManager::kDeltaPhiv0CrpH2,4,-1.*TMath::Pi(),TMath::Pi());

  // leg variables
  cf->AddVariable(AliDielectronVarManager::kPt,"0.8, 1.0, 1.1, 1.2, 1.5, 100.0",kTRUE);
  if(hasMC) cf->AddVariable(AliDielectronVarManager::kEta,"-0.9,0.9",kTRUE);
  //    cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,7,-1.5,5.5,kTRUE);
  //    cf->AddVariable(AliDielectronVarManager::kNclsITS,"1,2,3,4,5,6",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,"-3,-2.5,-2,2,2.5,3",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPio,"2.5,3.0,3.5,4.0,4.5,100",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kNclsTPC,"70, 90, 100, 120, 160",kTRUE);
  //    cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPro,"3.5,4.0,4.5,5.0,100",kTRUE);
  //    cf->AddVariable(AliDielectronVarManager::kTOFnSigmaEle,"-3,-2,2,3",kTRUE); break;
  //    cf->AddVariable(AliDielectronVarManager::kTRDpidQuality,"3.5, 4.5, 5.5, 6.5",kTRUE);
  //    if(!hasMC && isESD) cf->AddVariable(AliDielectronVarManager::kTRDchi2,"-1.,0.,2.,4.",kTRUE);
  */
  // mc steps
  if(hasMC) {
    if(cutDefinition==kTOFTRD) cf->SetStepForMCtruth();
    cf->SetStepsForMCtruthOnly();
    // cf->SetStepsForBackground();
  }

  die->SetCFManagerPair(cf);
}

void AddMCSignals(AliDielectron *die){
  //Do we have an MC handler?
  if (!hasMC) return;
  
  AliDielectronSignalMC* inclusiveJpsi = new AliDielectronSignalMC("inclusiveJpsi","Inclusive J/psi");
  inclusiveJpsi->SetLegPDGs(11,-11);
  inclusiveJpsi->SetMotherPDGs(443,443);
  inclusiveJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  inclusiveJpsi->SetFillPureMCStep(kTRUE);
  inclusiveJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  inclusiveJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(inclusiveJpsi);
  
  AliDielectronSignalMC* promptJpsi = new AliDielectronSignalMC("promptJpsi","Prompt J/psi");   // prompt J/psi (not from beauty decays)
  promptJpsi->SetLegPDGs(11,-11);
  promptJpsi->SetMotherPDGs(443,443);
  promptJpsi->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
  promptJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  promptJpsi->SetFillPureMCStep(kTRUE);
  promptJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  promptJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(promptJpsi);
  
  AliDielectronSignalMC* beautyJpsi = new AliDielectronSignalMC("beautyJpsi","Beauty J/psi");
  beautyJpsi->SetLegPDGs(11,-11);
  beautyJpsi->SetMotherPDGs(443,443);
  beautyJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  beautyJpsi->SetGrandMotherPDGs(500,500);
  beautyJpsi->SetFillPureMCStep(kTRUE);
  beautyJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  beautyJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  beautyJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(beautyJpsi);
  
  AliDielectronSignalMC* directJpsi = new AliDielectronSignalMC("directJpsi","Direct J/psi");   // embedded J/psi
  directJpsi->SetLegPDGs(11,-11);
  directJpsi->SetMotherPDGs(443,443);
  directJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  directJpsi->SetFillPureMCStep(kTRUE);
  directJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  directJpsi->SetMotherSources(AliDielectronSignalMC::kDirect, AliDielectronSignalMC::kDirect);
  directJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  directJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(directJpsi);
  
  AliDielectronSignalMC* conversionElePairs = new AliDielectronSignalMC("conversionElePairs","conversion electron pairs");  // pairs made from conversion (may be also from 2 different conversions)
  conversionElePairs->SetLegPDGs(11,-11);
  conversionElePairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  conversionElePairs->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  conversionElePairs->SetMotherPDGs(22,22);
  die->AddSignalMC(conversionElePairs);

  // prompt J/psi radiative channel
  AliDielectronSignalMC* promptJpsiRad = new AliDielectronSignalMC("promptJpsiRad","Prompt J/psi Radiative");   // prompt J/psi (not from beauty decays)
  promptJpsiRad->SetLegPDGs(11,-11);
  promptJpsiRad->SetMotherPDGs(443,443);
  promptJpsiRad->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
  promptJpsiRad->SetMothersRelation(AliDielectronSignalMC::kSame);
  promptJpsiRad->SetFillPureMCStep(kTRUE);
  promptJpsiRad->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  promptJpsiRad->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsiRad->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsiRad->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  promptJpsiRad->SetJpsiRadiative(AliDielectronSignalMC::kIsRadiative);
  //  die->AddSignalMC(promptJpsiRad);

  // prompt J/psi Non radiative channel
  AliDielectronSignalMC* promptJpsiNonRad = new AliDielectronSignalMC("promptJpsiNonRad","Prompt J/psi non-Radiative");   // prompt J/psi (not from beauty decays)
  promptJpsiNonRad->SetLegPDGs(11,-11);
  promptJpsiNonRad->SetMotherPDGs(443,443);
  promptJpsiNonRad->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
  promptJpsiNonRad->SetMothersRelation(AliDielectronSignalMC::kSame);
  promptJpsiNonRad->SetFillPureMCStep(kTRUE);
  promptJpsiNonRad->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  promptJpsiNonRad->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsiNonRad->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsiNonRad->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  promptJpsiNonRad->SetJpsiRadiative(AliDielectronSignalMC::kIsNotRadiative);
  //  die->AddSignalMC(promptJpsiNonRad);
}

void SetEtaCorrection()
{
  if (AliDielectronPID::GetEtaCorrFunction()) return;
  
  TString etaMap="$TRAIN_ROOT/jpsi_JPSI/EtaCorrMaps.root";
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if (trainRoot.IsNull()) etaMap="$ALICE_ROOT/PWGDQ/dielectron/files/EtaCorrMaps.root";
  if (gSystem->AccessPathName(gSystem->ExpandPathName(etaMap.Data()))){
    Error("ConfigPbPb","Eta map not found: %s",etaMap.Data());
    return;
  }

  TFile f(etaMap.Data());
  if (!f.IsOpen()) return;
  TList *keys=f.GetListOfKeys();

  for (Int_t i=0; i<keys->GetEntries(); ++i){
    TString kName=keys->At(i)->GetName();
    TPRegexp reg(kName);
    if (reg.MatchB(list)){
      printf(" Using Eta Correction Function: %s\n",kName.Data());
      AliDielectronPID::SetEtaCorrFunction((TF1*)f.Get(kName.Data()));
    }
  }
}

TVectorD *GetRunNumbers() {
  
  Double_t runLHC10h[] = { // all good runs based on RCT 29.Mai
    139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161, 137135
  };
  
  Double_t runLHC11h[] = { // all good runs based on RCT 29.Mai
    170593, 170572, 170388, 170387, 170315, 170313, 170312, 170311, 170309, 170308, 170306, 170270, 170269, 170268, 170230, 170228, 170207, 170204, 170203, 170193, 170163, 170159, 170155, 170091, 170089, 170088, 170085, 170084, 170083, 170081, 170040, 170027, 169965, 169923, 169859, 169858, 169855, 169846, 169838, 169837, 169835, 169591, 169590, 169588, 169587, 169586, 169557, 169555, 169554, 169553, 169550, 169515, 169512, 169506, 169504, 169498, 169475, 169420, 169419, 169418, 169417, 169415, 169411, 169238, 169167, 169160, 169156, 169148, 169145, 169144, 169138, 169099, 169094, 169091, 169045, 169044, 169040, 169035, 168992, 168988, 168826, 168777, 168514, 168512, 168511, 168467, 168464, 168460, 168458, 168362, 168361, 168342, 168341, 168325, 168322, 168311, 168310, 168115, 168108, 168107, 168105, 168076, 168069, 167988, 167987, 167985, 167920, 167915
  };
  
  // selection via environement variable (works only for gsi trains)

  
  if(list.Contains("LHC10h") || list.Contains("LHC11a10")) {
    Int_t size = (int) (sizeof(runLHC10h)/sizeof(Double_t));
    TVectorD *vec = new TVectorD(size+1);
    
    (*vec)[size] = runLHC10h[0] + 1;
    for (int i = 0; i < size; i++) {
      (*vec)[i] = runLHC10h[size-1-i];
    }
    //    vec->Print("");
    return vec;
  }

  if( list.Contains("LHC11h") || list.Contains("LHC12a17") ) {
    
    Int_t size = (int) (sizeof(runLHC11h)/sizeof(Double_t));
    TVectorD *vec = new TVectorD(size+1);
    
    (*vec)[size] = runLHC11h[0] + 1;
    for (int i = 0; i < size; i++) {
      (*vec)[i] = runLHC11h[size-1-i];
    }
    //   vec->Print("");
    return vec;
  }

  TVectorD *vec = new TVectorD(2);
  (*vec)[0] = 0;
  (*vec)[1] = 1;
  return vec;
     
}

TVectorD *GetDeltaPhiBins() {
  //
  // for in and out of event plane bins
  //
  Double_t pi = TMath::Pi();
  TVectorD *deltaPhi = new TVectorD(6);
  (*deltaPhi)[0] = -1.    *pi;
  (*deltaPhi)[1] = -3./4. *pi;
  (*deltaPhi)[2] = -1./4. *pi;
  (*deltaPhi)[3] = +1./4. *pi;
  (*deltaPhi)[4] = +3./4. *pi;
  (*deltaPhi)[5] = +1.    *pi;
  return deltaPhi;
}
