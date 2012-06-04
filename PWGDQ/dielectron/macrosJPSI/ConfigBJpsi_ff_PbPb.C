void InitHistograms(AliDielectron *die, Int_t cutDefinition);
void InitCF(AliDielectron* die, Int_t cutDefinition);

void SetupTrackCuts(Bool_t isESD, AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition);

void AddMCSignals(AliDielectron *die);
void SetEtaCorrection();
TVectorD *GetRunNumbers();

TString names=("TOFTRDany;TOFTRDfirst");
enum { kTOFTRD, kTOFTRD2};

TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntries();

Bool_t hasMC=kFALSE;

AliDielectron* ConfigBJpsi_ff_PbPb(Int_t cutDefinition, Bool_t isMC=kFALSE)
{
  //
  // Setup the instance of AliDielectron
  //
  

  // MC event handler?
  hasMC=isMC;
    //(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);    

  //ESD handler?
  Bool_t isESD=(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliESDInputHandler::Class());
  
  // switch off some configurations
  //switch(cutDefinition) {
   //case kTOFTRD:   
   //case kTOFTRD2:   
  //default:		return 0x0;      break;
  //}
  
  // create the actual framework object
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast()){
    name=arrNames->At(cutDefinition)->GetName();
  }
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()),
                                         Form("Track cuts: %s",name.Data()));
  // debugTree
  AliDielectronDebugTree *tree = new AliDielectronDebugTree(Form("%s",name.Data()),Form("Track cuts: %s",name.Data()));
  tree->SetOutputFileName(Form("jpsi_debug_tree_%s.root",name.Data())); 

  if(tree){
           tree->AddLegVariable(AliDielectronVarManager::kPt);
           tree->AddLegVariable(AliDielectronVarManager::kPx);
           tree->AddLegVariable(AliDielectronVarManager::kPy);
           tree->AddLegVariable(AliDielectronVarManager::kEta);
	   tree->AddLegVariable(AliDielectronVarManager::kXv);
           tree->AddLegVariable(AliDielectronVarManager::kYv);
           tree->AddLegVariable(AliDielectronVarManager::kE);
           tree->AddLegVariable(AliDielectronVarManager::kTPCnSigmaEle);
           tree->AddLegVariable(AliDielectronVarManager::kTPCsignal);
	   tree->AddLegVariable(AliDielectronVarManager::kTRDntracklets);
	   tree->AddLegVariable(AliDielectronVarManager::kNclsTRD);
	   tree->AddLegVariable(AliDielectronVarManager::kTOFnSigmaEle);
           tree->AddLegVariable(AliDielectronVarManager::kITSLayerFirstCls);
	   tree->AddPairVariable(AliDielectronVarManager::kM);
           tree->AddPairVariable(AliDielectronVarManager::kE);
           tree->AddPairVariable(AliDielectronVarManager::kP);
           tree->AddPairVariable(AliDielectronVarManager::kPt);
	   tree->AddPairVariable(AliDielectronVarManager::kY);
           tree->AddPairVariable(AliDielectronVarManager::kEta);
           tree->AddPairVariable(AliDielectronVarManager::kPairType);
	   tree->AddPairVariable(AliDielectronVarManager::kPseudoProperTime);
	   tree->AddPairVariable(AliDielectronVarManager::kPseudoProperTimeErr);
	   if(hasMC) tree->AddPairVariable(AliDielectronVarManager::kPseudoProperTimeResolution);
	   if(hasMC) tree->AddPairVariable(AliDielectronVarManager::kPseudoProperTimePull);
	   if(hasMC) tree->AddPairVariable(AliDielectronVarManager::kPdgCode);
           if(hasMC) tree->AddPairVariable(AliDielectronVarManager::kPdgCodeMother);
           if(hasMC) tree->AddPairVariable(AliDielectronVarManager::kPdgCodeGrandMother);
           if(hasMC) tree->AddPairVariable(AliDielectronVarManager::kIsJpsiPrimary);
	   tree->AddPairVariable(AliDielectronVarManager::kCentralitySPD);
	   tree->AddPairVariable(AliDielectronVarManager::kCentrality);
	   tree->AddPairVariable(AliDielectronVarManager::kNevents);
	   tree->AddPairVariable(AliDielectronVarManager::kXvPrim);
           tree->AddPairVariable(AliDielectronVarManager::kYvPrim); 
   	   tree->AddPairVariable(AliDielectronVarManager::kZvPrim);
           tree->AddPairVariable(AliDielectronVarManager::kXRes);
           tree->AddPairVariable(AliDielectronVarManager::kYRes);
	   tree->AddPairVariable(AliDielectronVarManager::kZRes);

	    }
  //
  // QA histogram setup
  //
  die->SetDebugTree(tree);


  // Monte Carlo Signals and TRD efficiency tables
  if(hasMC) {
    AddMCSignals(die);
    
    // trd tables
    TString pidTab="$TRAIN_ROOT/util/dielectron/dielectron/TRDpidEff_eleProb07_TRDntr4_6.root";
    TString trainRoot=gSystem->Getenv("TRAIN_ROOT");
    if (trainRoot.IsNull()) pidTab="$ALICE_ROOT/PWGDQ/dielectron/files/TRDpidEff_eleProb07_TRDntr4_6.root";

    if (gSystem->AccessPathName(gSystem->ExpandPathName(pidTab.Data())))
      Error("ConfigPbPb","PID table not found: %s",pidTab.Data());
    else 
      die->SetTRDcorrectionFilename(pidTab.Data());
  }
  
  // cut setup
  SetupTrackCuts(isESD,die,cutDefinition);
  SetupPairCuts(die,cutDefinition);
  
  // histogram setup
  if(cutDefinition == kTOFTRD  || 
	cutDefinition == kTOFTRD2   ) 
    InitHistograms(die,cutDefinition);
  
  // CF container setup
  //  InitCF(die,cutDefinition);
  
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
    
    
    }
    
  
  // prefilter settings
  //if(cutDefinition == kTOFTRD2) 
  //  die->SetPreFilterAllSigns();
  //else 
    die->SetPreFilterUnlikeOnly();
  
  // setup eta correction
  if(isESD) SetEtaCorrection();
  
  return die;
}

//______________________________________________________________________________________
void SetupTrackCuts(Bool_t isESD, AliDielectron *die, Int_t cutDefinition)
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
  
	// track cuts ESD and AOD
  AliDielectronVarCuts *varCuts = new AliDielectronVarCuts("VarCuts","VarCuts");
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.9,   0.9);
  //varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,     70.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.0);

   AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  if(isESD){
    varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
    switch(cutDefinition) {
    case kTOFTRD2: varCuts->AddCut(AliDielectronVarManager::kITSLayerFirstCls,-0.01,0.5); //ITS(0) = SPDfirst
      break;
    default:      varCuts->AddCut(AliDielectronVarManager::kITSLayerFirstCls,-0.01,1.5); //ITS(0-1) = SPDany
      break;
   }
  }else
      {
      switch(cutDefinition) {
      case kTOFTRD2: trkCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);  //ITS(0) = SPDfirst
      break;
      default:  trkCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);//ITS(0-1) = SPDany
      break;
      }
    }
 
  cuts->AddCut(varCuts);
 
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE);
  cuts->AddCut(trkCuts);
  
  //Do we have an MC handler?
  //  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronPID *pid = new AliDielectronPID("PID","PID");
  
  ////////////////////////////////// DATA
  if(!hasMC) {
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,3.5,0.,0.,kTRUE);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-100.,3.5,0.,0.,kTRUE);
    
    //if(cutDefinition==kTRD || cutDefinition>=kTOFTRD || cutDefinition>=kTOFTRD2) 
    if(cutDefinition>=kTOFTRD || cutDefinition>=kTOFTRD2) 
      pid->AddCut(AliDielectronPID::kTRDeleEff,AliPID::kElectron,.8,1.,3.5.,6.,kFALSE,
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
      
      TString list=gSystem->Getenv("LIST");
      
      //LHC11a10b
      if (list.Contains("LHC11a10b") || list.IsNull()) {
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
    if (list.Contains("LHC11a10b") || list.IsNull()) {
      nSigmaCorrection->SetPoint(0, 137161., -0.50-(0.28));
      nSigmaCorrection->SetPoint(1, 139510., -0.50-(0.28));
      pid->SetCorrGraph(nSigmaCorrection);
    }
    
  } //hasMC
  
  ////////////////////////////////// DATA + MC
  // pid cuts TPC + TOF & TRD
  pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.);
  if(cutDefinition==kTOFTRD || cutDefinition==kTOFTRD2) 
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3,3.,0.,0.,kFALSE,AliDielectronPID::kIfAvailable);
  
  //if(cutDefinition!=krec) 
  cuts->AddCut(pid);
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
  Double_t gCut = 0.05;             // default
  
  AliDielectronVarCuts *gammaCut=new AliDielectronVarCuts("gammaCut","gammaCut");
  gammaCut->AddCut(AliDielectronVarManager::kM,0.,gCut);
  die->GetPairPreFilter().AddCuts(gammaCut);
  
  
  // rapidity selection
  AliDielectronVarCuts *rapCut=new AliDielectronVarCuts("|Y|<.9","|Y|<.9");
  rapCut->AddCut(AliDielectronVarManager::kY,-0.9,0.9);
  die->GetPairFilter().AddCuts(rapCut);
  
}

//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //
  //  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(die->GetName(),die->GetTitle());
  
  
  //add histograms to event class
  histos->AddClass("Event");
  histos->UserHistogram("Event","VtxZ","Vertex Z;z (cm)",
                        300,-15.,15.,
                        AliDielectronVarManager::kZvPrim);
  //histos->UserHistogram("Event","Centrality","Centrality;centrality (%)",
  //                      "0.,5.,10.,20.,40.,50.,60.,80.,100.",
  //                      AliDielectronVarManager::kCentrality);
  histos->UserHistogram("Event","Centrality","Centrality;centrality (%)",
                        20,0.,100.,AliDielectronVarManager::kCentrality);

  histos->UserHistogram("Event","Multiplicity","Multiplicity V0;Multiplicity V0",
                        500,0.,25000., AliDielectronVarManager::kMultV0);
  histos->UserHistogram("Event","Cent_Mult","Centrality vs. Multiplicity;centrality (%);Multiplicity V0",
                        10,0.,100., 500,0.,25000.,AliDielectronVarManager::kCentrality,AliDielectronVarManager::kMultV0);
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
  
  if(cutDefinition <= kTOFTRD) {
    
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
  
   histos->UserHistogram("Pair","PseudoProperTime","Pseudoproper decay length; pseudoproper-decay-length[#mum];Entries/40#mum",
                          150,-0.3.,0.3,AliDielectronVarManager::kPseudoProperTime);


   }
  
  die->SetHistogramManager(histos);
}


void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setup the CF Manager if needed
  //
  //  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);  
  
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());
  
  // pair variables
  cf->AddVariable(AliDielectronVarManager::kM,125,0.,125*.04); //40Mev Steps
  if(cutDefinition!=kSubRndm) cf->AddVariable(AliDielectronVarManager::kPairType,11,0,11);
  
  if(cutDefinition <= kTOFTRD2) {
    
    // pair and event vars
      cf->AddVariable(AliDielectronVarManager::kCentrality,"0.,5.,10.,20.,40.,50.,60.,80.");
      cf->AddVariable(AliDielectronVarManager::kPt,"0., 1., 2.5, 5., 100.0");
      if(!hasMC) cf->AddVariable(AliDielectronVarManager::kZvPrim,20, -10., 10.);
      if(hasMC) cf->AddVariable(AliDielectronVarManager::kNacc,20,0.,3000.0);
      if(hasMC) cf->AddVariable(AliDielectronVarManager::kNVtxContrib,20,0.,4000.);
      if(hasMC) cf->AddVariable(AliDielectronVarManager::kTRDpidEffPair,101,0.0,1.01);
      //cf->AddVariable(AliDielectronVarManager::kY,"-0.8,0.8");
    
    //leg variables
      cf->AddVariable(AliDielectronVarManager::kPt,"0.8, 1.0, 1.1, 1.2, 1.5, 100.0",kTRUE);
      cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,"-3,-2.5,-2,2,2.5,3",kTRUE);
      //cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPio,"3.5,4.0,4.5,5.0,100",kTRUE);
      //    cf->AddVariable(AliDielectronVarManager::kTPCnSigmaPro,"3.5,4.0,4.5,5.0,100",kTRUE);
      
      // standard vars
        if(hasMC) cf->AddVariable(AliDielectronVarManager::kEta,"-0.9,0.9",kTRUE);
        //cf->AddVariable(AliDielectronVarManager::kNclsTPC,"70, 90, 100, 120, 160",kTRUE);
    }
    
  
  
  if(hasMC) {
    cf->AddVariable(AliDielectronVarManager::kRunNumber, GetRunNumbers() ); // LHC10h only
    //    cf->AddVariable(AliDielectronVarManager::kRunNumber, 170593-136831, 136831, 170593); // LHC10h -> LHC11h
    if(cutDefinition==kTOFTRD || cutDefinition==kTOFTRD2) cf->SetStepForMCtruth();
     cf->SetStepsForMCtruthOnly();  
	   // cf->SetStepsForBackground();   
  }
  
  die->SetCFManagerPair(cf);
}

void AddMCSignals(AliDielectron *die){
  //Do we have an MC handler?
  //Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
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
