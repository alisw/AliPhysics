void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition);

void AddMCSignals(AliDielectron *die);
void SetEtaCorrection();
TVectorD *GetRunNumbers();

AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition);

TString names=("TPC;TOF;TRD;rec;TOFTRD;TOFTRD2;ITScls;ITSamy;dca;chi2;Gam0;Gam01;Gam05;Gam10;Gam15;Gam20;EtaGap01;EtaGap02;EtaGap03;EtaGap04;EtaGap05;SubLS;SubRndm");
enum { kTPC=0, kTOF, kTRD, krec, kTOFTRD, kTOFTRD2, kITScls, kITSamy, kDCA, kChi, kGam0, kGam01, kGam05, kGam10, kGam15, kGam20, kEtaGap01, kEtaGap02, kEtaGap03, kEtaGap04, kEtaGap05, kSubLS, kSubRndm };

TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntries();

AliDielectron* ConfigJpsi_jb_PbPb(Int_t cutDefinition)
{
  //
  // Setup the instance of AliDielectron
  //
  
  // MC event handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);    

  //ESD handler?
  Bool_t isESD=(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliESDInputHandler::Class());
  
  // switch off some configurations
  switch(cutDefinition) {
    case kTPC:   
    case kTOF:   
    case kTRD:
      return 0x0;
      break;
      //case kTOFTRD:   
    case krec:   
      if(!hasMC) return 0x0;
      break;
      //case kTOFTRD2:   
    case kITScls:
    case kITSamy:
    case kDCA:
    case kChi:      
      return 0x0;
      break;
      //    case kGam0:   
      //    case kGam01:   
      //    case kGam05:   
      //    case kGam10:   
      //    case kGam15:   
      //    case kGam20:
    case kEtaGap01:   
    case kEtaGap02:   
    case kEtaGap03:   
    case kEtaGap04:   
    case kEtaGap05:   
    case kSubLS:
    case kSubRndm:
      if( hasMC) return 0x0;
      break;
  }
  
  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast()){
    name=arrNames->At(cutDefinition)->GetName();
  }
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()),
                                         Form("Track cuts: %s",name.Data()));
  
  // Monte Carlo Signals and TRD efficiency tables
  if(hasMC) {
    AddMCSignals(die);
    
    // trd tables
    TString pidTab="$TRAIN_ROOT/util/dielectron/dielectron/TRDpidEff_eleProb07_TRDntr4_6.root";
    TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
    if (trainRoot.IsNull()) pidTab="../PWGDQ/dielectron/files/TRDpidEff_eleProb07_TRDntr4_6.root";

    if (gSystem->AccessPathName(gSystem->ExpandPathName(pidTab.Data())))
      Error("ConfigPbPb","PID table not found: %s",pidTab.Data());
    else 
      die->SetTRDcorrectionFilename(pidTab.Data());
  }
  
  // cut setup
  SetupTrackCuts(die,cutDefinition);
  SetupPairCuts(die,cutDefinition);
  
  // histogram setup
  if(cutDefinition == kTOFTRD  || 
     cutDefinition == kGam0    || 
     cutDefinition == kTOFTRD2 || 
     cutDefinition >= kEtaGap01 ) 
    InitHistograms(die,cutDefinition);
  
  // CF container setup
  if(cutDefinition <  kEtaGap01 || 
     cutDefinition == kSubRndm  )
    InitCF(die,cutDefinition);
  
  // bgrd estimators
  if(!hasMC) {
    
    if(cutDefinition == kTOFTRD) {  
      // rotations
      AliDielectronTrackRotator *rot=new AliDielectronTrackRotator;
      rot->SetIterations(10);
      rot->SetConeAnglePhi(TMath::Pi());
      rot->SetStartAnglePhi(TMath::Pi());
      die->SetTrackRotator(rot);
      // mixing
      AliDielectronMixingHandler *mix=new AliDielectronMixingHandler;
      mix->AddVariable(AliDielectronVarManager::kZvPrim,20,-10.,10.);
      mix->AddVariable(AliDielectronVarManager::kCentrality,"0,5,10,20,50,80");
      mix->SetMixType(AliDielectronMixingHandler::kAll);
      mix->SetDepth(50);
      die->SetMixingHandler(mix);
    }
    
    
    // TPC event plane configurations
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
    poi->AddCut(AliDielectronVarManager::kM,2.92,3.16);     // particles of interest, jpsi mass window
    die->GetEventPlanePOIPreFilter().AddCuts(poi); 
    
    AliDielectronVarCuts *poiTrk = new AliDielectronVarCuts("PoItracks","PoItracks");
    poiTrk->AddCut(AliDielectronVarManager::kPt,0.8,1e30);
    die->GetEventPlanePreFilter().AddCuts(poiTrk);

    if(cutDefinition >= kEtaGap01 && 
       cutDefinition <  kSubLS     ) {
      AliDielectronVarCuts *etaGap = new AliDielectronVarCuts(AliDielectronVarManager::GetValueName(AliDielectronVarManager::kEta),"etaGap");
      etaGap->AddCut(AliDielectronVarManager::kEta,-1*gGap,gGap,kTRUE);
      die->GetEventPlanePreFilter().AddCuts(etaGap);
    }
    
    if(cutDefinition==kSubLS) die->SetLikeSignSubEvents();
    die->SetPreFilterEventPlane();
  }
  
  // prefilter settings
  if(cutDefinition == kTOFTRD2) 
    die->SetPreFilterAllSigns();
  else 
    die->SetPreFilterUnlikeOnly();
  
  // setup eta correction
  if(isESD) SetEtaCorrection();
  
  return die;
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
  
  //Pt cut, should make execution a bit faster
  AliDielectronVarCuts *pt = new AliDielectronVarCuts("Pt>.8","Pt>.8");
  pt->AddCut(AliDielectronVarManager::kPt,0.8,1e30);
  cuts->AddCut(pt);
  
  //ESD track cuts
  Bool_t isESD=(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliESDInputHandler::Class());
  if(isESD) cuts->AddCut(SetupESDtrackCuts(cutDefinition));
  
  // track cuts ESD and AOD
  AliDielectronVarCuts *varCuts = new AliDielectronVarCuts("VarCuts","VarCuts");
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.9,   0.9);
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,     70.0, 160.0);
  switch(cutDefinition) {
    case kTOFTRD2: varCuts->AddCut(AliDielectronVarManager::kITSLayerFirstCls,-0.01,0.5); //ITS(0) = SPDfirst
      break;
    default:       varCuts->AddCut(AliDielectronVarManager::kITSLayerFirstCls,-0.01,1.5); //ITS(0-1) = SPDany
      break;
  }
  cuts->AddCut(varCuts);
  
  AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE);
  cuts->AddCut(trkCuts);
  
  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronPID *pid = new AliDielectronPID("PID","PID");
  
  ////////////////////////////////// DATA
  if(!hasMC) {
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);   // merge with below ones vv
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,3.5,0.,0.,kTRUE);
    //if((cutDefinition!=kTOF && cutDefinition<kTOFTRD) || cutDefinition>=kGam0) 
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-100.,3.5,0.,0.,kTRUE);
    
    if(cutDefinition==kTRD || cutDefinition>=kTOFTRD || cutDefinition>=kTOFTRD2) 
      pid->AddCut(AliDielectronPID::kTRDeleEff,AliPID::kElectron,.8,1.,3.5.,6.,kFALSE,
                  AliDielectronPID::kIfAvailable,AliDielectronVarManager::kTRDpidQuality);
  }
  
  ////////////////////////////////// MC
  if(hasMC) {
    
    // electron
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);   //merge with data only ^^
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
      
      TString list=gSystem->Getenv("LIST");
      
      //LHC11a10b
      if (list.Contains("LHC11a10b")) {
        printf("LHC11a10b parameters\n");
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
  // pid cuts TOF & TRD
  if(cutDefinition==kTOF || cutDefinition>=kTOFTRD || cutDefinition>=kTOFTRD2) 
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3,3.,0.,0.,kFALSE,AliDielectronPID::kIfAvailable);
  
  if(cutDefinition!=krec) cuts->AddCut(pid);
  /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ PID CUTS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  
  
  // exclude conversion electrons selected by the tender
  AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
  noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
  cuts->AddCut(noconv);
  
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
  
  AliDielectronVarCuts *gammaCut=new AliDielectronVarCuts("gammaCut","gammaCut");
  gammaCut->AddCut(AliDielectronVarManager::kM,0.,gCut);
  die->GetPairPreFilter().AddCuts(gammaCut);
  
  
  // rapidity selection
  AliDielectronVarCuts *rapCut=new AliDielectronVarCuts("|Y|<.9","|Y|<.9");
  rapCut->AddCut(AliDielectronVarManager::kY,-0.9,0.9);
  die->GetPairFilter().AddCuts(rapCut);
  
}

//______________________________________________________________________________________
AliESDtrackCuts *SetupESDtrackCuts(Int_t cutDefinition)
{
  //
  // Setup default AliESDtrackCuts
  //
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;
  // reconstruction cuts
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);

  /*
   // basic track quality cuts  (basicQ)
   esdTrackCuts->SetMaxDCAToVertexZ(3.0);
   esdTrackCuts->SetMaxDCAToVertexXY(1.0);
   
   // acceptance cuts
   esdTrackCuts->SetEtaRange( -0.9 , 0.9 );
   //esdTrackCuts->SetPtRange( 0.8, 1e30 );

   // reconstruction cuts
   esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
   esdTrackCuts->SetRequireITSRefit(kTRUE);
   esdTrackCuts->SetRequireTPCRefit(kTRUE);
   
   esdTrackCuts->SetMinNClustersTPC(70);
   //  if(cutDefinition<=kChi || cutDefinition>=kEtaGap01) 
   //  else esdTrackCuts->SetMinNClustersTPC(120);
   esdTrackCuts->SetMaxChi2PerClusterTPC(4);
   //  if(cutDefinition!=kITSamy) esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);   //  SPD first
   // esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
   
   switch(cutDefinition) {
   case kTOFTRD2:     
   esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);  // SPDfirst
   break;
   default:
   esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);  // SPDany
   break;
   }
   */
  return esdTrackCuts;
}

//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(die->GetName(),die->GetTitle());
  
  
  //add histograms to event class
  histos->AddClass("Event");
  histos->UserHistogram("Event","VtxZ","Vertex Z;z (cm)",
                        300,-15.,15.,
                        AliDielectronVarManager::kZvPrim);
  histos->UserHistogram("Event","Centrality","Centrality;centrality (%)",
                        "0.,5.,10.,20.,40.,50.,60.,80.,100.",
                        AliDielectronVarManager::kCentrality);
  histos->UserHistogram("Event","Multiplicity","Multiplicity V0;Multiplicity V0",
                        500,0.,25000.,
                        AliDielectronVarManager::kMultV0);
  histos->UserHistogram("Event","Cent_Mult","Centrality vs. Multiplicity;centrality (%);Multiplicity V0",
                        10,0.,100., 500,0.,25000.,
                        AliDielectronVarManager::kCentrality,AliDielectronVarManager::kMultV0);
  histos->UserProfile("Event","Cent_Nacc",
                      "accepted tracks;centrality (%)",
                      AliDielectronVarManager::kNacc,
                      "0.,5.,10.,20.,40.,50.,60.,80.,100.",
                      AliDielectronVarManager::kCentrality);
  histos->UserProfile("Event","Cent_NVtxContrib",
                      "number of vertex contributors;centrality (%)",
                      AliDielectronVarManager::kNVtxContrib,
                      "0.,5.,10.,20.,40.,50.,60.,80.,100.",
                      AliDielectronVarManager::kCentrality);
  
  
  ////// FLOW //////
  if(0) {
  if(cutDefinition == kTOFTRD || cutDefinition >= kEtaGap01) {
    histos->UserHistogram("Event","TPCxH2","TPC Qx component;TPCxH2",
                          100,-1500.,1500.,
                          AliDielectronVarManager::kTPCxH2);
    histos->UserHistogram("Event","TPCyH2","TPC Qy component;TPCyH2",
                          100,-1500.,1500.,
                          AliDielectronVarManager::kTPCyH2);
    histos->UserHistogram("Event","TPCrpH2","TPC reaction plane; #Psi^{TPC}",
                          100,-2.,2.,
                          AliDielectronVarManager::kTPCrpH2);
    histos->UserHistogram("Event","TPCsub1xH2","TPC Qx component sub1;TPCsub1xH2",
                          100,-1500.,1500.,
                          AliDielectronVarManager::kTPCsub1xH2);
    histos->UserHistogram("Event","TPCsub1yH2","TPC Qy component sub1;TPCsub1yH2",
                          100,-1500.,1500.,
                          AliDielectronVarManager::kTPCsub1yH2);
    histos->UserHistogram("Event","TPCsub1rpH2","TPC reaction plane sub1; #Psi^{sub1}",
                          100,-2.,2.,
                          AliDielectronVarManager::kTPCsub1rpH2);
    histos->UserHistogram("Event","TPCsub2xH2","TPC Qx component sub2;TPCsub2xH2",
                          100,-1500.,1500.,
                          AliDielectronVarManager::kTPCsub2xH2);
    histos->UserHistogram("Event","TPCsub2yH2","TPC Qy component sub2;TPCsub2yH2",
                          100,-1500.,1500.,
                          AliDielectronVarManager::kTPCsub2yH2);
    histos->UserHistogram("Event","TPCsub2rpH2","TPC reaction plane sub2; #Psi^{sub2}",
                          100,-2.,2.,
                          AliDielectronVarManager::kTPCsub2rpH2);
    histos->UserHistogram("Event","TPCsub12DiffH2","TPC reaction plane diff; cos(2(#Psi^{sub1}-#Psi^{sub2}))",
                          100,-1.,1.,
                          AliDielectronVarManager::kTPCsub12DiffH2);
    /* // uncorrected eventplane
     histos->UserHistogram("Event","TPCxH2uc","TPC Qx component;TPCxH2uc",
     100,-1500.,1500.,
     AliDielectronVarManager::kTPCxH2uc);
     histos->UserHistogram("Event","TPCyH2uc","TPC Qy component;TPCyH2uc",
     100,-1500.,1500.,
     AliDielectronVarManager::kTPCyH2uc);
     histos->UserHistogram("Event","TPCrpH2uc","TPC reaction plane;TPCrpH2uc",
     100,-2.,2.,
     AliDielectronVarManager::kTPCrpH2uc);
     histos->UserHistogram("Event","TPCsub1xH2uc","TPC Qx component sub1;TPCsub1xH2uc",
     100,-1500.,1500.,
     AliDielectronVarManager::kTPCsub1xH2uc);
     histos->UserHistogram("Event","TPCsub1yH2uc","TPC Qy component sub1;TPCsub1yH2uc",
     100,-1500.,1500.,
     AliDielectronVarManager::kTPCsub1yH2uc);
     histos->UserHistogram("Event","TPCsub1rpH2uc","TPC reaction plane sub1;TPCsub1rpH2uc",
     100,-2.,2.,
     AliDielectronVarManager::kTPCsub1rpH2uc);
     histos->UserHistogram("Event","TPCsub2xH2uc","TPC Qx component sub2;TPCsub2xH2uc",
     100,-1500.,1500.,
     AliDielectronVarManager::kTPCsub2xH2uc);
     histos->UserHistogram("Event","TPCsub2yH2uc","TPC Qy component sub2;TPCsub2yH2uc",
     100,-1500.,1500.,
     AliDielectronVarManager::kTPCsub2yH2uc);
     histos->UserHistogram("Event","TPCsub2rpH2uc","TPC reaction plane sub2;TPCsub2rpH2uc",
     100,-2.,2.,
     AliDielectronVarManager::kTPCsub2rpH2uc);
     histos->UserHistogram("Event","TPCsub12DiffH2uc","TPC reaction plane difference;TPCsub12DiffH2uc",
     100,-1.,1.,
     AliDielectronVarManager::kTPCsub12DiffH2uc);
     */
    histos->UserHistogram("Event","V0ACrpH2","VZERO-AC RP; #Psi_{2}^{V0AC} (rad.)",
                          100,-2.0,2.0,
                          AliDielectronVarManager::kV0ACrpH2);
    histos->UserHistogram("Event","V0ArpH2","VZERO-A RP; #Psi_{2}^{V0A} (rad.)",
                          100,-2.0,2.0,
                          AliDielectronVarManager::kV0ArpH2);
    histos->UserHistogram("Event","V0CrpH2","VZERO-C RP; #Psi_{2}^{V0C} (rad.)",
                          100,-2.0,2.0,
                          AliDielectronVarManager::kV0CrpH2);
    
    histos->UserHistogram("Event","V0ATPCDiffH2","VZERO-A TPC diff; cos(2(#Psi^{V0A}-#Psi^{TPC}))",
                          300,-1.0,1.0,
                          AliDielectronVarManager::kV0ATPCDiffH2);
    histos->UserHistogram("Event","V0CTPCDiffH2","VZERO-C TPC diff; cos(2(#Psi^{V0C}-#Psi^{TPC}))",
                          300,-1.0,1.0,
                          AliDielectronVarManager::kV0CTPCDiffH2);
    histos->UserHistogram("Event","V0AV0CDiffH2","VZERO-A VZERO-C diff; cos(2(#Psi^{V0A}-#Psi^{V0C}))",
                          300,-1.0,1.0,
                          AliDielectronVarManager::kV0AV0CDiffH2);
    
    // centrality dependent event plane histograms
    histos->UserHistogram("Event","Cent_TPCrpH2","TPC RP;centrality (%);#Psi^{TPC} (rad.)",
                          10,0.,100.,100,-2.,2.,
                          AliDielectronVarManager::kCentrality,AliDielectronVarManager::kTPCrpH2);
    histos->UserHistogram("Event","Cent_TPCsub1rpH2","TPC-1 RP;centrality (%);#Psi^{sub1} (rad.)",
                          10,0.,100.,100,-2.,2.,
                          AliDielectronVarManager::kCentrality,AliDielectronVarManager::kTPCsub1rpH2);
    histos->UserHistogram("Event","Cent_TPCsub2rpH2","TPC-2 RP;centrality (%);#Psi^{sub2} (rad.)",
                          10,0.,100.,100,-2.,2.,
                          AliDielectronVarManager::kCentrality,AliDielectronVarManager::kTPCsub2rpH2);
    
    histos->UserHistogram("Event","Cent_V0ACrpH2","VZERO-AC RP;centrality (%);#Psi_{2}^{V0AC} (rad.)",
                          10,0.,100.,100,-2.0,2.0,
                          AliDielectronVarManager::kCentrality,AliDielectronVarManager::kV0ACrpH2);
    histos->UserHistogram("Event","Cent_V0ArpH2","VZERO-A RP;centrality (%);#Psi_{2}^{V0A} (rad.)",
                          10,0.,100.,100,-2.0,2.0,
                          AliDielectronVarManager::kCentrality,AliDielectronVarManager::kV0ArpH2);
    histos->UserHistogram("Event","Cent_V0CrpH2","VZERO-C RP;centrality (%);#Psi_{2}^{V0C} (rad.)",
                          10,0.,100.,100,-2.0,2.0,
                          AliDielectronVarManager::kCentrality,AliDielectronVarManager::kV0CrpH2);
    // for resolution calculation
    histos->UserHistogram("Event","Cent_V0ATPCDiffH2","VZERO-A TPC diff;centrality (%);cos(2(#Psi^{V0A}-#Psi^{TPC}))",
                          10,0.,100.,300,-1.0,1.0,
                          AliDielectronVarManager::kCentrality,AliDielectronVarManager::kV0ATPCDiffH2);
    histos->UserHistogram("Event","Cent_V0CTPCDiffH2","VZERO-C TPC diff;centrality (%);cos(2(#Psi^{V0C}-#Psi^{TPC}))",
                          10,0.,100.,300,-1.0,1.0,
                          AliDielectronVarManager::kCentrality,AliDielectronVarManager::kV0CTPCDiffH2);
    histos->UserHistogram("Event","Cent_V0AV0CDiffH2","VZERO-A VZERO-C diff;centrality (%);cos(2(#Psi^{V0A}-#Psi^{V0C}))",
                          10,0.,100.,300,-1.0,1.0,
                          AliDielectronVarManager::kCentrality,AliDielectronVarManager::kV0AV0CDiffH2);
    histos->UserHistogram("Event","Cent_TPCsub12DiffH2","TPC-sub1 TPC-sub2 diff;centrality (%);cos(2(#Psi^{sub1}-#Psi^{sub2}))",
                          10,0.,100.,300,-1.0,1.0,
                          AliDielectronVarManager::kCentrality,AliDielectronVarManager::kTPCsub12DiffH2);
    // detector effects
    histos->UserHistogram("Event","Cent_TPCsub12DiffH2Sin","TPC-sub1 TPC-sub2 diff;centrality (%);sin(2(#Psi^{sub1}-#Psi^{sub2}))",
                          10,0.,100.,300,-1.0,1.0,AliDielectronVarManager::kCentrality,AliDielectronVarManager::kTPCsub12DiffH2Sin);
    //// EPSelectionTask V0 information directly from the ESDs
    histos->UserHistogram("Event","Cent_v0ACrpH2","VZERO-AC RP;centrality (%);#Psi_{2}^{v0AC} (rad.)",
                          10,0.,100.,100,-2.0,2.0,
                          AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0ACrpH2);
    histos->UserHistogram("Event","Cent_v0ArpH2","VZERO-A RP;centrality (%);#Psi_{2}^{v0A} (rad.)",
                          10,0.,100.,100,-2.0,2.0,
                          AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0ArpH2);
    histos->UserHistogram("Event","Cent_v0CrpH2","VZERO-C RP;centrality (%);#Psi_{2}^{v0C} (rad.)",
                          10,0.,100.,100,-2.0,2.0,
                          AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0CrpH2);
    // for resolution calculation
    histos->UserHistogram("Event","Cent_v0ATPCDiffH2","VZERO-A TPC diff;centrality (%);cos(2(#Psi^{v0A}-#Psi^{TPC}))",
                          10,0.,100.,300,-1.0,1.0,
                          AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0ATPCDiffH2);
    histos->UserHistogram("Event","Cent_v0CTPCDiffH2","VZERO-C TPC diff;centrality (%);cos(2(#Psi^{v0C}-#Psi^{TPC}))",
                          10,0.,100.,300,-1.0,1.0,
                          AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0CTPCDiffH2);
    histos->UserHistogram("Event","Cent_v0Av0CDiffH2","VZERO-A VZERO-C diff;centrality (%);cos(2(#Psi^{v0A}-#Psi^{v0C}))",
                          10,0.,100.,300,-1.0,1.0,
                          AliDielectronVarManager::kCentrality,AliDielectronVarManager::kv0Av0CDiffH2);
    
  }
  }
  
  ////// MONTE CARLO //////
  /*
  if(cutDefinition == kTOFTRD && hasMC) {
    histos->AddClass("MCEvent");
    histos->UserHistogram("MCEvent","Cent_NJPsis","Centrality vs. generated incl. J/#psi per event;centrality (%);N_{J/#psi}",
                          10,0.,100., 21,-0.5,20.5,
                          AliDielectronVarManager::kCentrality,AliDielectronVarManager::kNumberOfJPsis);
  }
  */
  
  
  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");
  
  
  //Pair classes
  // to fill also mixed event histograms loop until 10
  for (Int_t i=0; i<3; ++i){
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
  }
  
  if(cutDefinition < kEtaGap01) {
    
    //legs from pair
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
    
    //add histograms to Track classes
    //histos->UserHistogram("Track","TOFbit","TOFbit;bit;#tracks",19,-9.5,9.5,AliDielectronVarManager::kTOFPIDBit);
    histos->UserHistogram("Track","Pt","Pt;Pt [GeV];#tracks",
                          400,0,20.,
                          AliDielectronVarManager::kPt);
    histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusteres;#tracks",
                          160,-0.5,159.5,
                          AliDielectronVarManager::kNclsTPC);
    histos->UserHistogram("Track","TPCsignalN","Number of Clusters TPC;TPC number clusteres;#tracks",
                          160,-0.5,159.5,
                          AliDielectronVarManager::kTPCsignalN);
    
    histos->UserHistogram("Track","dXY","dXY;dXY [cm];#tracks",
                          500,-1.,1.,
                          AliDielectronVarManager::kImpactParXY);
    histos->UserHistogram("Track","dZ","dZ;dZ [cm];#tracks",
                          600,-3.,3.,
                          AliDielectronVarManager::kImpactParZ);
    histos->UserHistogram("Track","Eta_Phi","Eta Phi Map; Eta; Phi;#tracks",
                          200,-1,1,200,0,6.285,
                          AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
    
    histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                          400,0.2,20.,200,0.,200.,
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);
    histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
                          400,0.2,20.,200,-10.,10.,
                          AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,kTRUE);
    
    histos->UserHistogram("Track","Ncl",";Number clusters TPC;Number clusters TPC",
                          160,-0.5,159.5,
                          AliDielectronVarManager::kNclsTPC);
    histos->UserHistogram("Track","NclFr",";Number of findable clusters (robust);Number findable clusters TPC",
                          160,-0.5,159.5,
                          AliDielectronVarManager::kNFclsTPCr);
    histos->UserHistogram("Track","Ncl_NclFr","Number of (findable) clusters TPC;found clusters;findable clusters",
                          160,-0.5,159.5,160,-0.5,159.5,
                          AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
    histos->UserHistogram("Track","NtrklTRD",";Number tracklets TRD for pid;Number tracklets TRD",
                          8,-0.5,7.5,
                          AliDielectronVarManager::kTRDpidQuality);
    
    //add histograms to Pair classes
    histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
                          300,.0,300*0.04,
                          AliDielectronVarManager::kM); // 40MeV bins, 12GeV/c2
    histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                          100,-1.,1.,
                          AliDielectronVarManager::kY);
    histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                          100,0.,3.15,
                          AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","Chi2NDF","#chi^{2}/NDF;#chi^{2}/NDF",
                          100,0.,20,
                          AliDielectronVarManager::kChi2NDF);
  }
  
  //// FLOW results use tprofiles
  if(0) {
  if(cutDefinition == kTOFTRD || cutDefinition == kTOFTRD2 || cutDefinition >= kEtaGap01) {
    
    histos->UserProfile("Pair","M_Cent_Pt_V0ACrpH2FlowV2",
                        "cos(2(#varphi-#Psi^{V0AC}));mass (GeV/c^{2});centrality (%);p_{T} (GeV/c)",
                        AliDielectronVarManager::kV0ACrpH2FlowV2,
                        125,0.,125*.04, 10, 0.,100., 200,0.,100.,
                        AliDielectronVarManager::kM, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kPt);
    
    histos->UserProfile("Pair","M_Cent_Pt_V0ArpH2FlowV2",
                        "cos(2(#varphi-#Psi^{V0A}));mass (GeV/c^{2});centrality (%);p_{T} (GeV/c)",
                        AliDielectronVarManager::kV0ArpH2FlowV2,
                        125,0.,125*.04, 10, 0.,100., 200,0.,100.,
                        AliDielectronVarManager::kM, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kPt);
    
    histos->UserProfile("Pair","M_Cent_Pt_V0CrpH2FlowV2",
                        "cos(2(#varphi-#Psi^{V0C}));mass (GeV/c^{2});centrality (%);p_{T} (GeV/c)",
                        AliDielectronVarManager::kV0CrpH2FlowV2,
                        125,0.,125*.04, 10, 0.,100., 200,0.,100.,
                        AliDielectronVarManager::kM, AliDielectronVarManager::kCentrality, AliDielectronVarManager::kPt);
    // 1D
    histos->UserProfile("Pair","M_V0ACrpH2FlowV2",
                        "cos(2(#varphi-#Psi^{V0AC}));mass (GeV/c^{2})",
                        AliDielectronVarManager::kV0ACrpH2FlowV2,
                        125,0.,125*.04,
                        AliDielectronVarManager::kM);
    
    histos->UserProfile("Pair","M_V0ArpH2FlowV2",
                        "cos(2(#varphi-#Psi^{V0AC}));mass (GeV/c^{2})",
                        AliDielectronVarManager::kV0ArpH2FlowV2,
                        125,0.,125*.04,
                        AliDielectronVarManager::kM);
    
    histos->UserProfile("Pair","M_V0CrpH2FlowV2",
                        "cos(2(#varphi-#Psi^{V0AC}));mass (GeV/c^{2})",
                        AliDielectronVarManager::kV0CrpH2FlowV2,
                        125,0.,125*.04,
                        AliDielectronVarManager::kM);
  }  
  }
  
  die->SetHistogramManager(histos);
}


void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setup the CF Manager if needed
  //
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);  
  
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());
  
  // pair variables
  cf->AddVariable(AliDielectronVarManager::kM,125,0.,125*.04); //40Mev Steps
  if(cutDefinition!=kSubRndm) cf->AddVariable(AliDielectronVarManager::kPairType,11,0,11);
  
  if(cutDefinition <  kGam0 || cutDefinition == kSubRndm) {
    
    // pair and event vars
    if(cutDefinition <= kChi || cutDefinition == kSubRndm) {
      cf->AddVariable(AliDielectronVarManager::kCentrality,"0.,5.,10.,20.,40.,50.,60.,80.");
      cf->AddVariable(AliDielectronVarManager::kPt,"0., 1., 2.5, 5., 100.0");
      if(!hasMC) cf->AddVariable(AliDielectronVarManager::kZvPrim,20, -10., 10.);
      if(hasMC) cf->AddVariable(AliDielectronVarManager::kNacc,20,0.,3000.0);
      if(hasMC) cf->AddVariable(AliDielectronVarManager::kNVtxContrib,20,0.,4000.);
      if(hasMC) cf->AddVariable(AliDielectronVarManager::kTRDpidEffPair,101,0.0,1.01);
      //cf->AddVariable(AliDielectronVarManager::kY,"-0.8,0.8");
    }
    
    //leg variables
    if(cutDefinition!=kSubRndm) {
      cf->AddVariable(AliDielectronVarManager::kPt,"0.8, 1.0, 1.1, 1.2, 1.5, 100.0",kTRUE);
      cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,"-3,-2.5,-2,2,2.5,3",kTRUE);
      //cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPio,"3.5,4.0,4.5,5.0,100",kTRUE);
      //    cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPro,"3.5,4.0,4.5,5.0,100",kTRUE);
      
      // standard vars
      if(cutDefinition<=kChi) {
        if(hasMC) cf->AddVariable(AliDielectronVarManager::kEta,"-0.9,0.9",kTRUE);
        //cf->AddVariable(AliDielectronVarManager::kNclsTPC,"70, 90, 100, 120, 160",kTRUE);
      }
    }
    
    switch(cutDefinition) {
      case kTPC:
      case krec:
      case kTOF: //cf->AddVariable(AliDielectronVarManager::kTOFnSigmaEle,"-3,-2,2,3",kTRUE); break;
      case kTRD: 
      case kTOFTRD: 
        // if(hasMC) cf->AddVariable(AliDielectronVarManager::kThetaCS,15,-1.,1.);
        //cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,7,-1.5,5.5,kTRUE); break;
        break;
      case kITScls: cf->AddVariable(AliDielectronVarManager::kNclsITS,"1,2,3,4,5,6",kTRUE);       break;
      case kITSamy: cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,7,-1.5,5.5,kTRUE); break;
      case kDCA:
        cf->AddVariable(AliDielectronVarManager::kITSLayerFirstCls,7,-1.5,5.5,kTRUE); 
        cf->AddVariable(AliDielectronVarManager::kImpactParXY,8,-2.,2.,kTRUE);
        cf->AddVariable(AliDielectronVarManager::kImpactParZ,8,-4.,4.,kTRUE); 
        break;
      case kChi:    cf->AddVariable(AliDielectronVarManager::kChi2NDF,"0,1,2,3,4,5",kTRUE); break;
        //    cf->AddVariable(AliDielectronVarManager::kTOFnSigmaEle,"-3,-2,2,3",kTRUE);
        //    cf->AddVariable(AliDielectronVarManager::kTOFPIDBit,"-.5,.5,1.5",kTRUE);
        //    cf->AddVariable(AliDielectronVarManager::kTRDpidQuality,"3.5, 4.5, 5.5, 6.5",kTRUE);
    }
    
  }
  
  
  if(hasMC) {
    cf->AddVariable(AliDielectronVarManager::kRunNumber, GetRunNumbers() ); // LHC10h -> LHC11h
    //    cf->AddVariable(AliDielectronVarManager::kRunNumber, 170593-136831, 136831, 170593); // LHC10h -> LHC11h
    if(cutDefinition==kTOFTRD || cutDefinition==kGam0) cf->SetStepForMCtruth();
    //    if(cutDefinition!=kTOFTRD) 
      cf->SetStepsForMCtruthOnly();  
    // cf->SetStepsForBackground();   
  }
  
  die->SetCFManagerPair(cf);
}

void AddMCSignals(AliDielectron *die){
  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
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
}

void SetEtaCorrection()
{
  if (AliDielectronPID::GetEtaCorrFunction()) return;
  
  TString list=gSystem->Getenv("LIST");

  TString etaMap="$TRAIN_ROOT/jpsi_JPSI/EtaCorrMaps.root";
  TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
  if (trainRoot.IsNull()) etaMap="../PWGDQ/dielectron/files/EtaCorrMaps.root";
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
      printf("Using Eta Correction Function: %s\n",kName.Data());
      AliDielectronPID::SetEtaCorrFunction((TF1*)f.Get(kName.Data()));
    }
  }
}

TVectorD *GetRunNumbers() {

  Double_t runLHC10h[] = { // all runs
    136851, 136854, 136879, 137042, 137045, 137124, 137125, 137132, 137133, 137135, 137136, 137137, 137161, 137162, 137163, 137165, 137230, 137231, 137232, 137235, 137236, 137243, 137365, 137366, 137370, 137430, 137431, 137432, 137434, 137439, 137440, 137441, 137443, 137530, 137531, 137539, 137541, 137544, 137546, 137549, 137595, 137608, 137609, 137638, 137639, 137685, 137686, 137689, 137691, 137692, 137693, 137704, 137718, 137722, 137724, 137748, 137751, 137752, 137843, 137844, 137847, 137848, 138125, 138126, 138150, 138151, 138153, 138154, 138190, 138192, 138197, 138200, 138201, 138225, 138275, 138359, 138364, 138396, 138438, 138439, 138442, 138469, 138533, 138534, 138578, 138579, 138582, 138583, 138620, 138621, 138624, 138637, 138638, 138652, 138653, 138662, 138666, 138730, 138731, 138732, 138736, 138737, 138740, 138742, 138795, 138796, 138826, 138828, 138830, 138831, 138836, 138837, 138870, 138871, 138872, 138924, 138965, 138972, 138973, 138976, 138977, 138978, 138979, 138980, 138982, 138983, 139024, 139025, 139028, 139029, 139030, 139031, 139034, 139036, 139037, 139038, 139042, 139104, 139105, 139107, 139110, 139172, 139173, 139308, 139309, 139310, 139311, 139314, 139316, 139328, 139329, 139360, 139437, 139438, 139439, 139440, 139441, 139465, 139466, 139467, 139470, 139471, 139503, 139504, 139505, 139507, 139510, 139511, 139513, 139514, 139517,
0
  };
  

  Int_t sizeLHC10h = (int) (sizeof(runLHC10h)/sizeof(Double_t)); 
  runLHC10h[sizeLHC10h-1] =   runLHC10h[sizeLHC10h-2] + 1.;
  TVectorD *vecLHC10h = new TVectorD(sizeLHC10h, runLHC10h);
  
  return vecLHC10h;

}
