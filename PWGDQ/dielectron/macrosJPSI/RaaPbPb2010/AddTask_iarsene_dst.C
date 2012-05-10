//__________________________________________________________________________________________
AliAnalysisTask *AddTask_iarsene_dst(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_iarsene_dst", "No analysis manager found.");
    return 0;
  }

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  //Get the current train configuration
  //TString trainConfig=gSystem->Getenv("CONFIG_FILE");
  
  gROOT->LoadMacro("AliCorrelationReducedEvent.cxx+g");
  gROOT->LoadMacro("AliAnalysisTaskCorrelationTree.cxx+g");
  
  //create task and add it to the manager
  AliAnalysisTaskCorrelationTree *task=new AliAnalysisTaskCorrelationTree("DSTTreeMaker");
  //if(trainConfig.Contains("PbPb")) task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral);
  //if(trainConfig=="pp") task->SetRejectPileup();
  if(!hasMC) task->UsePhysicsSelection(kFALSE);
  mgr->AddTask(task);
  
  
  task->SetEventFilter(CreateEventFilter());
  task->SetTrackFilter(CreateGlobalTrackFilter());
  task->SetFlowTrackFilter(CreateFlowTrackFilter());
  task->SetK0sPionCuts(CreateK0sPionCuts());
  task->SetLambdaProtonCuts(CreateLambdaProtonCuts());
  task->SetLambdaPionCuts(CreateLambdaPionCuts());
  task->SetK0sCuts(CreateK0sCuts());
  task->SetK0sMassRange(0.44,0.55);
  task->SetLambdaMassRange(1.090,1.14);
  task->SetLambdaCuts(CreateLambdaCuts());
  task->SetV0Histograms(CreateV0Histograms());
  
  task->AddDielectron(ConfigDielectron(0));   // J/psi -> e+e-
  task->AddDielectron(ConfigDielectron(1));   // phi -> K+K-
  
  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("diele_defaultTree",
                         TChain::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "diele_default");
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("qaHistos",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "dst_qaHistos.root");

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("dstTree",
                         TTree::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "dstTree.root");
  
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("friendTree",
                         TTree::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "dstTree_friend.root");
  cout << "output containers: " << endl
       << cOutputHist1 << endl
       << cOutputHist2 << endl
       << cOutputHist3 << endl;

  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  //mgr->ConnectOutput(task, 2, cOutputHist2);
  //mgr->ConnectOutput(task, 3, cOutputHist3);
  mgr->ConnectOutput(task, 2, cOutputHist3);
  
  return task;
}


//_______________________________________________________________________________________________
AliAnalysisCuts* CreateEventFilter() {
  //
  // Event wise cuts
  //
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  //eventCuts->SetRequireVertex();
  //eventCuts->SetMinVtxContributors(1);
  //eventCuts->SetVertexZ(-10.,10.);
  return eventCuts;
}


//_______________________________________________________________________________________________
AliAnalysisCuts* CreateGlobalTrackFilter() {
  //
  // Cuts for tracks to be written in the dst
  //
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  // general ESD cuts ---------------------------------
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;
  // basic track quality cuts  (basicQ)
  //esdTrackCuts->SetMaxDCAToVertexZ(10.0);
  //esdTrackCuts->SetMaxDCAToVertexXY(3.0);
  //esdTrackCuts->SetEtaRange( -0.9 , 0.9 );
  //esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  //  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetPRange(.1,1e30);
  //esdTrackCuts->SetMinNClustersTPC(60);
  //esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  //  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  //  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  cuts->AddCut(esdTrackCuts);
  
  return cuts;
}


//_______________________________________________________________________________________________
AliAnalysisCuts* CreateFlowTrackFilter() {
  //
  // Cuts for tracks to be used for the event plane q-vector
  // These cuts are applied aditionally to the global track cuts
  //
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;
  esdTrackCuts->SetPtRange(0.2,2.0);
  esdTrackCuts->SetEtaRange(-0.8, 0.8);
  esdTrackCuts->SetMinNClustersTPC(70);
  return esdTrackCuts;
}


//_______________________________________________________________________________________________
AliAnalysisCuts* CreateK0sPionCuts() {
  //
  // Cuts on the K0s pions (tracking cuts, pid cuts, ...) 
  //
  AliESDtrackCuts *pionCuts = new AliESDtrackCuts;
  pionCuts->SetPtRange(0.15,100.0);
  return pionCuts;
}


//_______________________________________________________________________________________________
AliAnalysisCuts* CreateLambdaPionCuts() {
  //
  // Cuts on the Lambda pions (tracking cuts, pid cuts, ...) 
  //
  AliESDtrackCuts *pionCuts = new AliESDtrackCuts;
  pionCuts->SetPtRange(0.15,100.0);
  return pionCuts;
}


//_______________________________________________________________________________________________
AliAnalysisCuts* CreateLambdaProtonCuts() {
  //
  // Cuts on the Lambda protons (tracking cuts, pid cuts, ...) 
  //
  AliESDtrackCuts *protonCuts = new AliESDtrackCuts;
  protonCuts->SetPtRange(0.15,100.0);
  return protonCuts;
}


//______________________________________________________________________________________
AliESDv0Cuts* CreateK0sCuts() {
  //
  // Cuts on the V0s with K0s hypothesis
  //
  
  AliESDv0Cuts *cuts=new AliESDv0Cuts();
  //cuts->SetMinDcaPosToVertex(-1);
  //cuts->SetMinDcaNegToVertex(-1);
  //cuts->SetMaxChi2(10);
  //cuts->SetMaxDcaV0Daughters(0.3);
  //cuts->SetMinRadius(3.0);
  //cuts->SetMaxRadius(90.0);
  //cuts->SetMinCosinePointingAngle(0.9);
  //cuts->SetRequireOnFlyStatus(kTRUE);
  //cuts->SetMaxDcaV0ToVertex(0.5);
  //cuts->SetPRange(0.,1.0e10);
  //cuts->SetPtRange(0.,1.0e10);
  return cuts;
}


//______________________________________________________________________________________
AliESDv0Cuts* CreateLambdaCuts() {
  //
  // Cuts on the V0s with Lambda hypothesis
  //
  
  AliESDv0Cuts *cuts=new AliESDv0Cuts();
  //cuts->SetMinDcaPosToVertex(-1);
  //cuts->SetMinDcaNegToVertex(-1);
  //cuts->SetMaxChi2(10);
  //cuts->SetMaxDcaV0Daughters(0.3);
  //cuts->SetMinRadius(3.0);
  //cuts->SetMaxRadius(90.0);
  //cuts->SetMinCosinePointingAngle(0.9);
  //cuts->SetRequireOnFlyStatus(kTRUE);
  //cuts->SetMaxDcaV0ToVertex(0.5);
  //cuts->SetPRange(0.,1.0e10);
  //cuts->SetPtRange(0.,1.0e10);
  return cuts;
}


//______________________________________________________________________________________
AliDielectronHistos* CreateV0Histograms()
{
  //
  // Initialise the histograms
  //
  
  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos("V0Histograms","");
        
  //add histograms to Track classes ---------------------------------------------------------------------------
  histos->SetReservedWords("V0Track");
  
  //Track classes
  histos->AddClass("V0Track_Pos"); histos->AddClass("V0Track_Neg"); 
  
  // kinematic acceptance histograms
  histos->UserHistogram("V0Track","Pt","Pt;Pt (GeV/c);#tracks",400,0,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("V0Track","Eta_P","Eta_P;#eta;P (GeV/c)",40, -1.0, +1.0, 200,0,10., AliDielectronVarManager::kEta, AliDielectronVarManager::kP);
  histos->UserHistogram("V0Track","Eta_TRDPhi","Eta Phi Map; Eta; Phi;#tracks",
                        40,-1,1,200,-3.15,3.15,AliDielectronVarManager::kEta,AliDielectronVarManager::kTRDphi);  

  // DCA diagnostic histograms
  histos->UserHistogram("V0Track","dXY","dXY;dXY (cm);#tracks",500,-1.,1.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("V0Track","dZ","dZ;dZ (cm);#tracks",600,-3.,3.,AliDielectronVarManager::kImpactParZ);
  histos->UserHistogram("V0Track","dZ_dXY","dZ_dXY;dZ (cm);dXY (cm)",600,-3.,3.,500,-1.,1.,AliDielectronVarManager::kImpactParZ, AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("V0Track","P_dXY","P_dXY;P (GeV/c);dXY (cm)",200,0.,10.,500,-1.,1.,AliDielectronVarManager::kP, AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("V0Track","P_dZ","P_dZ;P (GeV/c);dZ (cm)",200,0.,10.,600,-3.,3.,AliDielectronVarManager::kP, AliDielectronVarManager::kImpactParZ);

  // ITS
  histos->UserHistogram("V0Track","ITSsignal","ITS signal;ITS dE/dx",1000,0.0,1000.0, AliDielectronVarManager::kITSsignal);  
  histos->UserHistogram("V0Track","ITSsignal_P","ITS signal vs P;P (GeV/c); ITS dE/dx",300, 0.0, 3.0, 200,0.0,1000.0, AliDielectronVarManager::kP, AliDielectronVarManager::kITSsignal);
  
  // TPC
  histos->UserHistogram("V0Track","TPCchi2","#chi^{2}/cluster;#chi^{2}/cluster;#tracks",200,0.0,5.0, AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("V0Track","Pt_TPCchi2","#chi^{2}/cluster vs pt;p_{T} (GeV/c);#chi^{2}/cluster",400,0.,20., 200,0.0,5.0,
			AliDielectronVarManager::kPt,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("V0Track","TPCnCls","Number of Clusters TPC;TPC number clusters;#tracks",160,-0.5,159.5, AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("V0Track","Phi_TPCnCls","Number of Clusters TPC vs #phi;#phi (rad.);TPC number clusters",200,0,6.285, 160,-0.5,159.5,
			AliDielectronVarManager::kPhi,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("V0Track","Pin_TPCnCls","Number of Clusters TPC vs TPC inner param P;p (GeV/c);TPC number clusters",400,0,20., 160,-0.5,159.5,
			AliDielectronVarManager::kPIn, AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("V0Track","TPCnCls_dEdx","Number of Clusters TPC vs dE/dx; TPC number clusters; dE/dx (arb.units)", 160,-0.5,159.5, 200,0,200.,
			AliDielectronVarManager::kNclsTPC, AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("V0Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                        400,0.2,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("V0Track","dEdx_Eta","dEdx;#eta;TPC signal (arb units);#tracks",
                        110,-1.1,1.1,200,0.,200.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("V0Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
                        400,0.2,20.,200,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("V0Track","TPCnSigmaEle_Eta","TPC number of sigmas Electrons;#eta;TPC number of sigmas Electrons;#tracks",
                        110,-1.1,1.1,200,-5.,5.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);  
  histos->UserHistogram("V0Track","TPCnSigmaPio_P","TPC number of sigmas Pions;P [GeV];TPC number of sigmas Pions;#tracks",
			400,0.2,20.,200,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
  histos->UserHistogram("V0Track","TPCnSigmaPio_Eta","TPC number of sigmas Pions;#eta;TPC number of sigmas Pions;#tracks",
                        110,-1.1,1.1,200,-5.,5.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaPio);  
  histos->UserHistogram("V0Track","TPCnSigmaKao_P","TPC number of sigmas Kaons;P [GeV];TPC number of sigmas Kaons;#tracks",
			400,0.2,20.,200,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
  histos->UserHistogram("V0Track","TPCnSigmaKao_Eta","TPC number of sigmas Kaons;#eta;TPC number of sigmas Kaons;#tracks",
                        110,-1.1,1.1,200,-5.,5.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaKao);  
  histos->UserHistogram("V0Track","TPCnSigmaPro_P","TPC number of sigmas Protons;P [GeV];TPC number of sigmas Protons;#tracks",
			400,0.2,20.,200,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);
  histos->UserHistogram("V0Track","TPCnSigmaPro_Eta","TPC number of sigmas Protons;#eta;TPC number of sigmas Protons;#tracks",
                        110,-1.1,1.1,200,-5.,5.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaPro);  
  histos->UserHistogram("V0Track","POut","Track POut;P (GeV/c)",100,0.0,10.0, AliDielectronVarManager::kPOut);

  // TOF
  histos->UserHistogram("V0Track","Pin_TOFsignal","TOF signal vs TPC inner param P; P [GeV]; TOF signal",
			120,0.0,6.,1200,0.,120000.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFsignal);
  histos->UserHistogram("V0Track","Pin_TOFbeta","TOF #beta vs TPC inner param P; P [GeV]; TOF #beta",
			120,0.0,6.,120,0.,1.2,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
  histos->UserHistogram("V0Track","TOFnSigmaPio_Pin","TOF number of sigmas Pions;P [GeV];TOF number of sigmas Pions;#tracks",
			400,0.2,20.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPio);
  histos->UserHistogram("V0Track","TOFnSigmaPio_eta","TOF number of sigmas Pions vs #eta;#eta;TOF number of sigmas Pions;#tracks",
			100,-1.2,1.2,200,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaPio);
  histos->UserHistogram("V0Track","TOFnSigmaKao_Pin","TOF number of sigmas Kaons;P [GeV];TOF number of sigmas Kaons;#tracks",
			400,0.2,20.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao);
  histos->UserHistogram("V0Track","TOFnSigmaKao_eta","TOF number of sigmas Kaons vs #eta;#eta;TOF number of sigmas Kaons;#tracks",
			100,-1.2,1.2,200,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaKao);
  histos->UserHistogram("V0Track","TOFnSigmaPro_Pin","TOF number of sigmas Protons;P [GeV];TOF number of sigmas Protons;#tracks",
			400,0.2,20.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro);
  histos->UserHistogram("V0Track","TOFnSigmaPro_eta","TOF number of sigmas Protons vs #eta;#eta;TOF number of sigmas Protons;#tracks",
			100,-1.2,1.2,200,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaPro);
  return histos;
}


//____________________________________________________________________________________________
AliDielectron* ConfigDielectron(Int_t cutDefinition)
{
  //
  // Setup the instance of an AliDielectron
  //
  const Char_t* dieleNames[2] = {"JpsiToEE", "PhiToKK"};
    
  AliDielectron *die =
    new AliDielectron(dieleNames[cutDefinition],
                      Form("Track cuts: %s", dieleNames));
  //die->SetEstimatorFilename("$TRAIN_ROOT/util/dielectron/dielectron/estimators.root");
  if(cutDefinition==0) {       // J/psi -> e+e-
    die->SetLegPdg(11,11);
    die->SetMotherPdg(443);
  }
  if(cutDefinition==1) {       // phi -> K+K-
    die->SetLegPdg(321,321);
    die->SetMotherPdg(333);
  }
  
  // cut setup
  SetupDielectronTrackCuts(die,cutDefinition);
  SetupDielectronPairCuts(die,cutDefinition);

  // Set MC signals
  SetDielectronMCSignals(die);
  
  //
  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielectron framework histograms will be filled
  //
  InitDielectronHistograms(die,cutDefinition);
 
  //setup eta correction
  SetEtaCorrection();

  return die;
}

//______________________________________________________________________________________
void SetupDielectronTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //

  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  die->GetTrackFilter().AddCuts(cuts);
  //ESD quality cuts  
  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;
  // basic track quality cuts  (basicQ)
  esdTrackCuts->SetMaxDCAToVertexZ(3.0);
  esdTrackCuts->SetMaxDCAToVertexXY(1.0);
  if(cutDefinition==1) {
    esdTrackCuts->SetMaxDCAToVertexZ(0.3);
    esdTrackCuts->SetMaxDCAToVertexXY(0.3);
  }
  if(cutDefinition==1) esdTrackCuts->SetEtaRange(-1.0,1.0);
  if(cutDefinition==0) esdTrackCuts->SetEtaRange(-0.9,0.9);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  //esdTrackCuts->SetRequireITSRefit(kTRUE);
  if(cutDefinition==0) esdTrackCuts->SetRequireTPCRefit(kTRUE);
  if(cutDefinition==0) esdTrackCuts->SetPRange(.9,1e30);   // for J/psi
  if(cutDefinition==1) esdTrackCuts->SetPRange(.1,1e30);   // for phi
  if(cutDefinition==0) esdTrackCuts->SetMinNClustersTPC(60);
  if(cutDefinition==0) esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  //esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  cuts->AddCut(esdTrackCuts);
    
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  AliDielectronPID *pid = new AliDielectronPID("PID","PID cut");
  if(cutDefinition==0) {   // J/psi->ee, phi->ee
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.0, 4.5, 0.0, 0.0,  kFALSE, AliDielectronPID::kRequire); // TPC 3-sigma inclusion for electron
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,  -100.0, 2.0, 0.0, 10.0, kTRUE,  AliDielectronPID::kRequire); // TPC 3.0-sigma exclusion for proton
//    pid->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3.0, 3.0, 0.0, 10.0, kTRUE,  AliDielectronPID::kRequire); // TPC 3.0-sigma exclusion for kaon
//    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -3.0, 3.0, 0.0, 10.0, kTRUE,  AliDielectronPID::kRequire); // TPC 3.0-sigma exclusion for pion
  }
  if(cutDefinition==1) {     // phi -> KK
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kKaon,   -3.5, 3.5, 0.6, 100.0, kFALSE, AliDielectronPID::kRequire, AliDielectronVarManager::kPIn);  // TOF 3-sigma inclusion for kaon
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,   -3.0, 3.5, 0.3, 100.0, kFALSE, AliDielectronPID::kRequire, AliDielectronVarManager::kPIn);  // TPC 3-sigma inclusion for kaon
    /*    pid->AddCut(AliDielectronPID::kTPC,AliPID::kProton, -2.0, 2.0, 0.0, 0.0, kTRUE,  AliDielectronPID::kRequire);  // TPC 2.0-sigma exclusion for  proton
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,   -2.0, 2.0, 0.0, 0.0, kTRUE,  AliDielectronPID::kRequire);  // TPC 2.0-sigma exclusion for  pion
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kProton, -2.0, 2.0, 0.6, 3.0, kTRUE,  AliDielectronPID::kIfAvailable);  // TOF 2.0-sigma exclusion for  proton
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kPion,   -2.0, 2.0, 0.6, 3.0, kTRUE,  AliDielectronPID::kIfAvailable);  // TOF 2.0-sigma exclusion for  pion
*/
  }
    
  // ++++++++++++++++++++++++++++++++++++++++
  // shifts for the nSigma electrons
  TGraph* nSigmaCorrection = new TGraph();
  // b period (LHC10f7a)
  nSigmaCorrection->SetPoint(0, 114000., -0.23-(0.0));
  nSigmaCorrection->SetPoint(1, 117222., -0.23-(0.0));
  // c period (LHC10f7a)
  nSigmaCorrection->SetPoint(2, 119000., -0.23-(0.0));
  nSigmaCorrection->SetPoint(3, 120829., -0.23-(0.0));
  // d period (LHC10f7a)
  nSigmaCorrection->SetPoint(4, 121000., -0.353-(-0.111));
  nSigmaCorrection->SetPoint(5, 127000., -0.353-(-0.111));
  // e period (LHC10f7a)
  nSigmaCorrection->SetPoint(6, 127500., -0.351-(-0.116));
  nSigmaCorrection->SetPoint(7, 130900., -0.351-(-0.116));
  // LHC11a period (LHC11b10b)
  nSigmaCorrection->SetPoint(8, 146680., -0.70-(+0.381));   //
  nSigmaCorrection->SetPoint(9, 146870., -0.70-(+0.381));   //
  if(hasMC)
    pid->SetCorrGraph(nSigmaCorrection);
  
  cuts->AddCut(pid);
}

//______________________________________________________________________________________
void SetupDielectronPairCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the pair cuts
  //
  
  //Invariant mass selection
  AliDielectronVarCuts *cuts=new AliDielectronVarCuts("|y|<0.9","|Y|<.9");
  if(cutDefinition==0) cuts->AddCut(AliDielectronVarManager::kY,-0.9,0.9);
  if(cutDefinition==1) cuts->AddCut(AliDielectronVarManager::kY,-1.0,1.0);
  if(cutDefinition==1) cuts->AddCut(AliDielectronVarManager::kM, 0.95,1.1);
  die->GetPairFilter().AddCuts(cuts);

  // Configure prefilter to remove conversions
  AliDielectronVarCuts *gammaCut=new AliDielectronVarCuts("gammaCut","gammaCut");
  gammaCut->AddCut(AliDielectronVarManager::kM,0.,.05);
  if(cutDefinition==0) {
    die->GetPairPreFilter().AddCuts(gammaCut);
    die->SetPreFilterUnlikeOnly();
  }  
}

//______________________________________________________________________________________
void InitDielectronHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //
  
  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(die->GetName(),die->GetTitle());
  
  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");
  
  //Track classes
  //to fill also track info from 2nd event loop until 2
  for (Int_t i=0; i<2; ++i){
    histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
  }
  
  //Pair classes
  // to fill also mixed event histograms loop until 10
  for (Int_t i=0; i<3; ++i){
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
  }
  
  //legs from pair
  for (Int_t i=0; i<3; ++i){
    histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i)));
  }
    
  //add histograms to event class
  if (cutDefinition==0) {
    histos->AddClass("Event");
    histos->UserHistogram("Event","VtxZ","Vertex Z;Z[cm]",500,-25.,25.,AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","NTrk","NTrk",500,-0.5,499.5,AliDielectronVarManager::kNTrk);
    histos->UserHistogram("Event","VtxZ_NaccTrcklts","VtxZ_NaccTrcklts;VtxZ;NaccTrcklts", 
			  240,-12.,12.,200,-0.5,199.5,AliDielectronVarManager::kZvPrim, AliDielectronVarManager::kNaccTrcklts);
 
    histos->UserHistogram("Event","MultV0A","VZERO mult",2000,0.,2000.,AliDielectronVarManager::kMultV0A);
    histos->UserHistogram("Event","MultV0C","VZERO mult",2000,0.,2000.,AliDielectronVarManager::kMultV0C);
    histos->UserHistogram("Event","MultV0","VZERO mult",2000,0.,2000.,AliDielectronVarManager::kMultV0);
    histos->UserHistogram("Event","AdcV0A","VZERO adc",10000,0.,100000.,AliDielectronVarManager::kAdcV0A);
    histos->UserHistogram("Event","AdcV0C","VZERO adc",10000,0.,100000.,AliDielectronVarManager::kAdcV0C);
    histos->UserHistogram("Event","AdcV0","VZERO adc",10000,0.,100000.,AliDielectronVarManager::kAdcV0);
    histos->UserHistogram("Event","MultV0A_MultV0C","MultV0A_MultV0C;MultV0A;MultV0C", 
			  2000,0.0,2000.,2000,0.0,2000,AliDielectronVarManager::kMultV0A, AliDielectronVarManager::kMultV0C);
    histos->UserHistogram("Event","AdcV0A_AdcV0C","AdcV0A_AdcV0C;AdcV0A;AdcV0C", 
			  1000,0.0,20000.,1000,0.0,20000,AliDielectronVarManager::kAdcV0A, AliDielectronVarManager::kAdcV0C);
    histos->UserHistogram("Event","MultV0A_AdcV0A","MultV0A_AdcV0A;MultV0A;AdcV0A", 
			  2000,0.0,2000.,2000,0.0,20000,AliDielectronVarManager::kMultV0A, AliDielectronVarManager::kAdcV0A);
    histos->UserHistogram("Event","MultV0C_AdcV0C","MultV0C_AdcV0C;MultV0C;AdcV0C", 
			  2000,0.0,2000.,2000,0.0,20000,AliDielectronVarManager::kMultV0C, AliDielectronVarManager::kAdcV0C);
  }
    
  //add histograms to Track classes ---------------------------------------------------------------------------
  // kinematic acceptance histograms
  histos->UserHistogram("Track","Pt","Pt;Pt (GeV/c);#tracks",400,0,20.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","Eta_P","Eta_P;#eta;P (GeV/c)",40, -1.0, +1.0, 200,0,10., AliDielectronVarManager::kEta, AliDielectronVarManager::kP);
  histos->UserHistogram("Track","Eta_Phi","Eta Phi Map; Eta; Phi;#tracks",
                        40,-1,1,200,0,6.285,AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);  

  // DCA diagnostic histograms
  histos->UserHistogram("Track","dXY","dXY;dXY (cm);#tracks",500,-1.,1.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ","dZ;dZ (cm);#tracks",600,-3.,3.,AliDielectronVarManager::kImpactParZ);
  histos->UserHistogram("Track","dZ_dXY","dZ_dXY;dZ (cm);dXY (cm)",600,-3.,3.,500,-1.,1.,AliDielectronVarManager::kImpactParZ, AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","P_dXY","P_dXY;P (GeV/c);dXY (cm)",200,0.,10.,500,-1.,1.,AliDielectronVarManager::kP, AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","P_dZ","P_dZ;P (GeV/c);dZ (cm)",200,0.,10.,600,-3.,3.,AliDielectronVarManager::kP, AliDielectronVarManager::kImpactParZ);

    // ITS
  histos->UserHistogram("Track","ITSsignal","ITS signal;ITS dE/dx",1000,0.0,1000.0, AliDielectronVarManager::kITSsignal);  
  histos->UserHistogram("Track","ITSsignal_P","ITS signal vs P;P (GeV/c); ITS dE/dx",300, 0.0, 3.0, 400,0.0,1000.0, AliDielectronVarManager::kP, AliDielectronVarManager::kITSsignal);
  
  // TPC
  histos->UserHistogram("Track","TPCchi2","#chi^{2}/cluster;#chi^{2}/cluster;#tracks",200,0.0,5.0, AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","Pt_TPCchi2","#chi^{2}/cluster vs pt;p_{T} (GeV/c);#chi^{2}/cluster",400,0.,20., 200,0.0,5.0,
			AliDielectronVarManager::kPt,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","TPCnCls","Number of Clusters TPC;TPC number clusters;#tracks",160,-0.5,159.5, AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","Phi_TPCnCls","Number of Clusters TPC vs #phi;#phi (rad.);TPC number clusters",200,0,6.285, 160,-0.5,159.5,
			AliDielectronVarManager::kPhi,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","Pin_TPCnCls","Number of Clusters TPC vs TPC inner param P;p (GeV/c);TPC number clusters",400,0,20., 160,-0.5,159.5,
			AliDielectronVarManager::kPIn, AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","TPCnCls_dEdx","Number of Clusters TPC vs dE/dx; TPC number clusters; dE/dx (arb.units)", 160,-0.5,159.5, 200,0,200.,
			AliDielectronVarManager::kNclsTPC, AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","dEdx_P","dEdx;P [GeV];TPC signal (arb units);#tracks",
                        400,0.2,20.,200,0.,200.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","TPCnSigmaEle_P","TPC number of sigmas Electrons;P [GeV];TPC number of sigmas Electrons;#tracks",
                        400,0.2,20.,200,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaPio_P","TPC number of sigmas Pions;P [GeV];TPC number of sigmas Pions;#tracks",
			400,0.2,20.,200,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
  histos->UserHistogram("Track","TPCnSigmaKao_P","TPC number of sigmas Kaons;P [GeV];TPC number of sigmas Kaons;#tracks",
			400,0.2,20.,200,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
  histos->UserHistogram("Track","TPCnSigmaPro_P","TPC number of sigmas Protons;P [GeV];TPC number of sigmas Protons;#tracks",
			400,0.2,20.,200,-5.,5.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);
  histos->UserHistogram("Track","dEdx_Eta","dEdx;#eta;TPC signal (arb units);#tracks",
                        100,-1.0,1.0,200,0.,200.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","TPCnSigmaEle_Eta","TPC number of sigmas Electrons;#eta;TPC number of sigmas Electrons;#tracks",
                        100,-1.0,1.0,200,-5.,5.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaPio_Eta","TPC number of sigmas Pions;#eta;TPC number of sigmas Pions;#tracks",
			100,-1.0,1.0,200,-5.,5.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaPio);
  histos->UserHistogram("Track","TPCnSigmaKao_Eta","TPC number of sigmas Kaons;#eta;TPC number of sigmas Kaons;#tracks",
			100,-1.0,1.0,200,-5.,5.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaKao);
  histos->UserHistogram("Track","TPCnSigmaPro_Eta","TPC number of sigmas Protons;#eta;TPC number of sigmas Protons;#tracks",
			100,-1.0,1.0,200,-5.,5.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaPro);

  // TOF
  histos->UserHistogram("Track","Pin_TOFsignal","TOF signal vs TPC inner param P; P [GeV]; TOF signal",
			120,0.0,6.,1200,0.,120000.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFsignal);
  histos->UserHistogram("Track","Pin_TOFbeta","TOF #beta vs TPC inner param P; P [GeV]; TOF #beta",
			120,0.0,6.,120,0.,1.2,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
  histos->UserHistogram("Track","TOFnSigmaPio_Pin","TOF number of sigmas Pions;P [GeV];TOF number of sigmas Pions;#tracks",
			400,0.2,20.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPio);
  histos->UserHistogram("Track","TOFnSigmaPio_eta","TOF number of sigmas Pions vs #eta;#eta;TOF number of sigmas Pions;#tracks",
			100,-1.2,1.2,200,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaPio);
  histos->UserHistogram("Track","TOFnSigmaKao_Pin","TOF number of sigmas Kaons;P [GeV];TOF number of sigmas Kaons;#tracks",
			400,0.2,20.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao);
  histos->UserHistogram("Track","TOFnSigmaKao_eta","TOF number of sigmas Kaons vs #eta;#eta;TOF number of sigmas Kaons;#tracks",
			100,-1.2,1.2,200,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaKao);
  histos->UserHistogram("Track","TOFnSigmaPro_Pin","TOF number of sigmas Protons;P [GeV];TOF number of sigmas Protons;#tracks",
			400,0.2,20.,200,-10.,10.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro);
  histos->UserHistogram("Track","TOFnSigmaPro_eta","TOF number of sigmas Protons vs #eta;#eta;TOF number of sigmas Protons;#tracks",
			100,-1.2,1.2,200,-10.,10.,AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaPro);

  //add histograms to Pair classes
  Int_t nMassBins = 10000;
  Double_t massRange[2] = {-0.005, 100.005};
  if(cutDefinition==1) {   // Phi -> K+K-
    nMassBins = 4000;
    massRange[0] = 0.8; massRange[1] = 1.2;
  }
  histos->UserHistogram("Pair","InvMass","Inv.Mass;Inv. Mass [GeV];#pairs",
			nMassBins, massRange[0], massRange[1], AliDielectronVarManager::kM);
  histos->UserHistogram("Pair","InvMass_Pt","Inv.Mass vs pt;Inv. Mass [GeV]; pt (GeV/c);#pairs",
			nMassBins/10, massRange[0], massRange[1], 100, 0.0, 50.0, AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","InvMass_Cent","Inv.Mass vs centrality;Inv. Mass [GeV]; centrality;#pairs",
			nMassBins/10, massRange[0], massRange[1], 20, 0.0, 100.0, AliDielectronVarManager::kM, AliDielectronVarManager::kCentrality);
  histos->UserHistogram("Pair","InvMass_Rap","Inv.Mass vs rapidity;Inv. Mass [GeV]; rapidity;#pairs",
			nMassBins/10, massRange[0], massRange[1], 100, -1.0, 1.0, AliDielectronVarManager::kM, AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","InvMass_OpAngle","Inv.Mass vs opening angle;Inv. Mass [GeV]; op. angle (rad.);#pairs",
			nMassBins/10, massRange[0], massRange[1], 100, 0.0, 4.0, AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","InvMass_ThetaHE","Inv.Mass vs cos #theta^{*}_{HE};Inv. Mass [GeV]; cos #theta^{*}_{HE};#pairs",
			nMassBins/10, massRange[0], massRange[1], 100, -1.0, 1.0, AliDielectronVarManager::kM, AliDielectronVarManager::kThetaHE);
  histos->UserHistogram("Pair","InvMass_ThetaCS","Inv.Mass vs cos #theta^{*}_{CS};Inv. Mass [GeV]; cos #theta^{*}_{CS};#pairs",
			nMassBins/10, massRange[0], massRange[1], 100, -1.0, 1.0, AliDielectronVarManager::kM, AliDielectronVarManager::kThetaCS);
  histos->UserHistogram("Pair","Pt","Transverse momentum; P_{t} [GeV/c];#pairs",
			1000, 0.0, 50.0, AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","Pt_Y","Transverse momentum vs rapidity; rapidity;  P_{t} [GeV/c];#pairs",
			100, -1.0, 1.0, 200, 0.0, 50.0, AliDielectronVarManager::kY, AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","Pt_ThetaCS","Transverse momentum vs cos #theta^{*}_{CS}; cos #theta^{*}_{CS};  P_{t} [GeV/c];#pairs",
			100, -1.0, 1.0, 200, 0.0, 50.0, AliDielectronVarManager::kThetaCS, AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","Pt_ThetaHE","Transverse momentum vs cos #theta^{*}_{HE}; cos #theta^{*}_{HE};  P_{t} [GeV/c];#pairs",
			100, -1.0, 1.0, 200, 0.0, 50.0, AliDielectronVarManager::kThetaHE, AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","Pt_Phi","Transverse momentum vs phi; #phi (rad.);  P_{t} [GeV/c];#pairs",
			200, -6.3, 6.3, 200, 0.0, 50.0, AliDielectronVarManager::kPhi, AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair","Phi","#phi; #phi (rad.);#pairs",
			1000, -6.3, 6.3, AliDielectronVarManager::kPhi);
  histos->UserHistogram("Pair","Rapidity","Rapidity;Rapidity;#pairs",
                        100,-1.,1.,AliDielectronVarManager::kY);
  histos->UserHistogram("Pair","Rapidity_OpAngle","Rapidity vs opening angle;Rapidity; Op.Angle (rad.);#pairs",
                        100,-1.,1.,200, 0.0, 4.0, AliDielectronVarManager::kY, AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","Rapidity_ThetaCS","Rapidity vs cos #theta^{*}_{CS};Rapidity; cos #theta^{*}_{CS};#pairs",
                        100,-1.,1.,200, -1.0, 1.0, AliDielectronVarManager::kY, AliDielectronVarManager::kThetaCS);
  histos->UserHistogram("Pair","Rapidity_ThetaHE","Rapidity vs cos #theta^{*}_{HE};Rapidity; cos #theta^{*}_{HE};#pairs",
                        100,-1.,1.,200, -1.0, 1.0, AliDielectronVarManager::kY, AliDielectronVarManager::kThetaHE);
  histos->UserHistogram("Pair","OpeningAngle","Opening angle;angle",
                        100,0.,3.15,AliDielectronVarManager::kOpeningAngle);
  histos->UserHistogram("Pair","Radius","Radius;R[cm]",
                        1000,0.,300.,AliDielectronVarManager::kR);
  histos->UserHistogram("Pair","ThetaHE","cos #theta^{*}_{HE}; cos #theta^{*}_{HE}",
                        220,-1.1,1.1,AliDielectronVarManager::kThetaHE);
  histos->UserHistogram("Pair","PhiHE","#varphi^{*}_{HE};#varphi^{*}_{HE} (rad.)",
                        160,-3.2,3.2,AliDielectronVarManager::kPhiHE);
  histos->UserHistogram("Pair","ThetaCS","cos #theta^{*}_{CS}; cos #theta^{*}_{CS}",
                        220,-1.1,1.1,AliDielectronVarManager::kThetaCS);
  histos->UserHistogram("Pair","PhiCS","#varphi^{*}_{CS};#varphi^{*}_{CS} (rad.)",
                        160,-3.2,3.2,AliDielectronVarManager::kPhiCS);
  histos->UserHistogram("Pair","Lxy","Pseudo-proper time;Lxy (cm.)",
                        1000,0.0,2.0,AliDielectronVarManager::kPseudoProperTime);
  histos->UserHistogram("Pair","Chi2NDF","Pair #chi^{2}/NDF; #chi^{2}/NDF",
                        1000,0.0,10.0,AliDielectronVarManager::kChi2NDF);
  
  
  die->SetHistogramManager(histos);
}


//_____________________________________________________________________________________________
void SetDielectronMCSignals(AliDielectron *die)
{
  // J/psi sources (direct + feed down from charm higher states)
  // 0
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

  // 1
  AliDielectronSignalMC* beautyJpsi = new AliDielectronSignalMC("beautyJpsi","beauty hadron -> J/psi");  // J/psi->e+e- from beauty hadron decays
  beautyJpsi->SetLegPDGs(11,-11);
  beautyJpsi->SetMotherPDGs(443,443);
  beautyJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  beautyJpsi->SetGrandMotherPDGs(503,503);
  beautyJpsi->SetFillPureMCStep(kTRUE);
  beautyJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  beautyJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  beautyJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  beautyJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(beautyJpsi);

  // 2
  AliDielectronSignalMC* beautyMesonJpsi = new AliDielectronSignalMC("beautyMesonJpsi","beauty meson -> J/psi");  // J/psi->e+e- from beauty meson decays
  beautyMesonJpsi->SetLegPDGs(11,-11);
  beautyMesonJpsi->SetMotherPDGs(443,443);
  beautyMesonJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  beautyMesonJpsi->SetGrandMotherPDGs(500,500);
  beautyMesonJpsi->SetFillPureMCStep(kTRUE);
  beautyMesonJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  beautyMesonJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  beautyMesonJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  beautyMesonJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(beautyMesonJpsi);
  
  // 3
  AliDielectronSignalMC* openBeautyMesonJpsi = new AliDielectronSignalMC("openBeautyMesonJpsi","open beauty meson -> J/psi");  // J/psi->e+e- from open beauty meson decays
  openBeautyMesonJpsi->SetLegPDGs(11,-11);
  openBeautyMesonJpsi->SetMotherPDGs(443,443);
  openBeautyMesonJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  openBeautyMesonJpsi->SetGrandMotherPDGs(501,501);
  openBeautyMesonJpsi->SetFillPureMCStep(kTRUE);
  openBeautyMesonJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  openBeautyMesonJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  openBeautyMesonJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  openBeautyMesonJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(openBeautyMesonJpsi);

  // 4
  AliDielectronSignalMC* chic0Jpsi = new AliDielectronSignalMC("chic0Jpsi","chic0 -> J/psi + X");  // J/psi->e+e- from chic0 decays
  chic0Jpsi->SetLegPDGs(11,-11);
  chic0Jpsi->SetMotherPDGs(443,443);
  chic0Jpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  chic0Jpsi->SetGrandMotherPDGs(10441,10441);
  chic0Jpsi->SetFillPureMCStep(kTRUE);
  chic0Jpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  chic0Jpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  chic0Jpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  chic0Jpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(chic0Jpsi);
  
  // 5
  AliDielectronSignalMC* chic1Jpsi = new AliDielectronSignalMC("chic1Jpsi","chic1 -> J/psi + X");  // J/psi->e+e- from chic1 decays
  chic1Jpsi->SetLegPDGs(11,-11);
  chic1Jpsi->SetMotherPDGs(443,443);
  chic1Jpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  chic1Jpsi->SetGrandMotherPDGs(20443,20443);
  chic1Jpsi->SetFillPureMCStep(kTRUE);
  chic1Jpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  chic1Jpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  chic1Jpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  chic1Jpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(chic1Jpsi);
  
  // 6
  AliDielectronSignalMC* chic2Jpsi = new AliDielectronSignalMC("chic2Jpsi","chic2 -> J/psi + X");  // J/psi->e+e- from chic2 decays
  chic2Jpsi->SetLegPDGs(11,-11);
  chic2Jpsi->SetMotherPDGs(443,443);
  chic2Jpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  chic2Jpsi->SetGrandMotherPDGs(445,445);
  chic2Jpsi->SetFillPureMCStep(kTRUE);
  chic2Jpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  chic2Jpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  chic2Jpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  chic2Jpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(chic2Jpsi);
  
  // 7
  AliDielectronSignalMC* psiPrimeJpsi = new AliDielectronSignalMC("psiPrimeJpsi","psi(2S) -> J/psi + X");  // J/psi->e+e- from psi(2S) decays
  psiPrimeJpsi->SetLegPDGs(11,-11);
  psiPrimeJpsi->SetMotherPDGs(443,443);
  psiPrimeJpsi->SetMothersRelation(AliDielectronSignalMC::kSame);
  psiPrimeJpsi->SetGrandMotherPDGs(445,445);
  psiPrimeJpsi->SetFillPureMCStep(kTRUE);
  psiPrimeJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  psiPrimeJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  psiPrimeJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  psiPrimeJpsi->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(psiPrimeJpsi);

  // physical backgrounds (electrons pairs from other sources and their combinatorics)
  // 8
  AliDielectronSignalMC* diEleContinuum = new AliDielectronSignalMC("diEleContinuum","di-electron continuum");     // all di-electrons originating in the collision
  diEleContinuum->SetLegPDGs(11,-11);
  diEleContinuum->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleContinuum->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleContinuum->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleContinuum);

  // 9
  AliDielectronSignalMC* diEleCharm = new AliDielectronSignalMC("diEleCharm","di-electrons from charm");  // dielectrons originating from charm hadrons (not neccessary from same mother)
  diEleCharm->SetLegPDGs(11,-11);
  diEleCharm->SetMotherPDGs(403,403);
  diEleCharm->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleCharm->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleCharm->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleCharm->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleCharm);

  // 10
  AliDielectronSignalMC* diEleOpenCharm = new AliDielectronSignalMC("diEleOpenCharm","di-electrons from open charm");  // dielectrons originating from open charm hadrons
  diEleOpenCharm->SetLegPDGs(11,-11);
  diEleOpenCharm->SetMotherPDGs(402,402);
  diEleOpenCharm->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleOpenCharm->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenCharm->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOpenCharm->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleOpenCharm);

  // 11
  AliDielectronSignalMC* diEleOpenCharmJpsi = new AliDielectronSignalMC("diEleOpenCharmJpsi","1 leg from open charm + 1 jpsi leg");  // 1 leg from open charm + 1 leg from jpsi
  diEleOpenCharmJpsi->SetLegPDGs(11,-11);
  diEleOpenCharmJpsi->SetMotherPDGs(402,443);
  diEleOpenCharmJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleOpenCharmJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenCharmJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOpenCharmJpsi->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleOpenCharmJpsi);
  
  // 12
  AliDielectronSignalMC* diEleBeauty = new AliDielectronSignalMC("diEleBeauty","di-electrons from beauty");  // dielectrons originating from beauty hadrons (not neccessary from same mother)
  diEleBeauty->SetLegPDGs(11,-11);
  diEleBeauty->SetMotherPDGs(503,503);
  diEleBeauty->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleBeauty->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleBeauty->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleBeauty->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleBeauty);

  // 13
  AliDielectronSignalMC* diEleOpenBeauty = new AliDielectronSignalMC("diEleOpenBeauty","di-electrons from open beauty");  // dielectrons originating from open beauty hadrons
  diEleOpenBeauty->SetLegPDGs(11,-11);
  diEleOpenBeauty->SetMotherPDGs(502,502);
  diEleOpenBeauty->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleOpenBeauty->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenBeauty->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOpenBeauty->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleOpenBeauty);

  // 14
  AliDielectronSignalMC* diEleBeautyJpsi = new AliDielectronSignalMC("diEleBeautyJpsi","1 leg from beauty + 1 jpsi leg");  // 1 leg from beauty + 1 leg from jpsi
  diEleBeautyJpsi->SetLegPDGs(11,-11);
  diEleBeautyJpsi->SetMotherPDGs(503,443);
  diEleBeautyJpsi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleBeautyJpsi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleBeautyJpsi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleBeautyJpsi->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleBeautyJpsi);

  // 15
  AliDielectronSignalMC* diEleBeautyOpenCharm = new AliDielectronSignalMC("diEleBeautyOpenCharm","1 leg from beauty + 1 leg from open charm");  // 1 leg from open charm + 1 leg from beauty
  diEleBeautyOpenCharm->SetLegPDGs(11,-11);
  diEleBeautyOpenCharm->SetMotherPDGs(503,403);
  diEleBeautyOpenCharm->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleBeautyOpenCharm->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleBeautyOpenCharm->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleBeautyOpenCharm->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleBeautyOpenCharm);
  
  // 16
  AliDielectronSignalMC* diEleStrange = new AliDielectronSignalMC("diEleStrange","di-electrons from strange particles");  // dielectrons originating from strange hadrons (not neccessary from same mother)
  diEleStrange->SetLegPDGs(11,-11);
  diEleStrange->SetMotherPDGs(300,300);
  diEleStrange->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleStrange->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleStrange->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleStrange->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleStrange);
  
  // 17
  AliDielectronSignalMC* diEleStrangeCharm = new AliDielectronSignalMC("diEleStrangeCharm","1 leg from strange + 1 leg from charm");  // 1 leg from strange + 1 leg from charm
  diEleStrangeCharm->SetLegPDGs(11,-11);
  diEleStrangeCharm->SetMotherPDGs(300,403);
  diEleStrangeCharm->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleStrangeCharm->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleStrangeCharm->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleStrangeCharm->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleStrangeCharm);
  
  // 18
  AliDielectronSignalMC* diEleStrangeBeauty = new AliDielectronSignalMC("diEleStrangeBeauty","1 leg from strange + 1 leg from beauty");  // 1 leg from strange + 1 leg from beauty
  diEleStrangeBeauty->SetLegPDGs(11,-11);
  diEleStrangeBeauty->SetMotherPDGs(300,503);
  diEleStrangeBeauty->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleStrangeBeauty->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleStrangeBeauty->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleStrangeBeauty->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleStrangeBeauty);
  
  // 19
  AliDielectronSignalMC* diEleLight1 = new AliDielectronSignalMC("diEleLight1","1 leg from light particles + 1 leg from everything");  // 1 leg from light hadrons (100+ PYTHIA codes) + 1 leg from everything
  diEleLight1->SetLegPDGs(11,-11);
  diEleLight1->SetMotherPDGs(100,0);
  diEleLight1->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleLight1->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleLight1->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleLight1->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleLight1);
  
  // 20
  AliDielectronSignalMC* diEleLight1Charm = new AliDielectronSignalMC("diEleLight1Charm","1 leg from light particles + 1 leg from charm");  // 1 leg from light hadrons (100+ PYTHIA codes) + 1 leg from charm
  diEleLight1Charm->SetLegPDGs(11,-11);
  diEleLight1Charm->SetMotherPDGs(100,403);
  diEleLight1Charm->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleLight1Charm->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleLight1Charm->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleLight1Charm->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleLight1Charm);
  
  // 21
  AliDielectronSignalMC* diEleLight1Beauty = new AliDielectronSignalMC("diEleLight1Beauty","1 leg from light particles + 1 leg from beauty");  // 1 leg from light hadrons (100+ PYTHIA codes) + 1 leg from beauty
  diEleLight1Beauty->SetLegPDGs(11,-11);
  diEleLight1Beauty->SetMotherPDGs(100,503);
  diEleLight1Beauty->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleLight1Beauty->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleLight1Beauty->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleLight1Beauty->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleLight1Beauty);
  
  // 22
  AliDielectronSignalMC* diEleLight2 = new AliDielectronSignalMC("diEleLight2","1 leg from light particles + 1 leg from everything");  // 1 leg from light hadrons (200+ PYTHIA codes) + 1 leg from everything
  diEleLight2->SetLegPDGs(11,-11);
  diEleLight2->SetMotherPDGs(200,0);
  diEleLight2->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleLight2->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleLight2->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleLight2->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleLight2);
  
  // 23
  AliDielectronSignalMC* diEleLight2Charm = new AliDielectronSignalMC("diEleLight2Charm","1 leg from light particles + 1 leg from charm");  // 1 leg from light hadrons (200+ PYTHIA codes) + 1 leg from charm
  diEleLight2Charm->SetLegPDGs(11,-11);
  diEleLight2Charm->SetMotherPDGs(100,403);
  diEleLight2Charm->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleLight2Charm->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleLight2Charm->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleLight2Charm->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleLight2Charm);
  
  // 24
  AliDielectronSignalMC* diEleLight2Beauty = new AliDielectronSignalMC("diEleLight2Beauty","1 leg from light particles + 1 leg from beauty");  // 1 leg from light hadrons (100+ PYTHIA codes) + 1 leg from beauty
  diEleLight2Beauty->SetLegPDGs(11,-11);
  diEleLight2Beauty->SetMotherPDGs(100,503);
  diEleLight2Beauty->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  diEleLight2Beauty->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleLight2Beauty->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleLight2Beauty->SetFillPureMCStep(kTRUE);
  die->AddSignalMC(diEleLight2Beauty);
  
  // background from secondary electrons
  // 25
  AliDielectronSignalMC* secondaryElectrons = new AliDielectronSignalMC("secondaryElectrons","Secondary electrons");   // all di-electrons from secondary electrons (interaction with detector)
  secondaryElectrons->SetLegPDGs(11,-11);
  secondaryElectrons->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  secondaryElectrons->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(secondaryElectrons);

  // 26
  AliDielectronSignalMC* primarySecElePairs = new AliDielectronSignalMC("primarySecElePairs","Primary+Secondary electron pairs");  // primary-secondary pairs
  primarySecElePairs->SetLegPDGs(11,-11);
  primarySecElePairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  primarySecElePairs->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kSecondary);
  die->AddSignalMC(primarySecElePairs);

  // 27
  AliDielectronSignalMC* conversionElePairs = new AliDielectronSignalMC("conversionElePairs","conversion electron pairs");      // pairs made from conversion (may be also from 2 different conversions)
  conversionElePairs->SetLegPDGs(11,-11);
  conversionElePairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  conversionElePairs->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  conversionElePairs->SetMotherPDGs(22,22);
  die->AddSignalMC(conversionElePairs);

  // misidentification
  // 28
  AliDielectronSignalMC* allEleMisIdPairs = new AliDielectronSignalMC("allEleMisIdPairs","all electron+misid. pairs");  // one true electron + a mis-id electron (all sources included)
  allEleMisIdPairs->SetLegPDGs(11,11,kFALSE,kTRUE);
  allEleMisIdPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(allEleMisIdPairs);

  // 29
  AliDielectronSignalMC* allMisIdMisIdPairs = new AliDielectronSignalMC("allMisIdMisIdPairs","all misid.+misid. pairs");  // mis-id + mis-id
  allMisIdMisIdPairs->SetLegPDGs(11,11,kTRUE,kTRUE);
  allMisIdMisIdPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(allMisIdMisIdPairs);

  // 30
  AliDielectronSignalMC* elePionPairs = new AliDielectronSignalMC("elePionPairs","electron+pion pairs");    // true electron + mis-id pion
  elePionPairs->SetLegPDGs(11,211);
  elePionPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(elePionPairs);

  // 31
  AliDielectronSignalMC* eleProtonPairs = new AliDielectronSignalMC("eleProtonPairs","Electron+proton pairs");  // true electron + mis-id proton
  eleProtonPairs->SetLegPDGs(11,2212);
  eleProtonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(eleProtonPairs);
  
  // 32
  AliDielectronSignalMC* eleKaonPairs = new AliDielectronSignalMC("eleKaonPairs","electron+kaon pairs");   // true electron + mis-id kaon
  eleKaonPairs->SetLegPDGs(11,321);
  eleKaonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(eleKaonPairs);

  // 33
  AliDielectronSignalMC* eleProtonPairs = new AliDielectronSignalMC("eleProtonPairs","Electron+proton pairs");  // true electron + mis-id proton
  eleProtonPairs->SetLegPDGs(11,2212);
  eleProtonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(eleProtonPairs);

  // 34
  AliDielectronSignalMC* piPiPairs = new AliDielectronSignalMC("piPiPairs","pion+pion pairs");    // mis-id pion + mis-id pion
  piPiPairs->SetLegPDGs(211,211);
  piPiPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(piPiPairs);

  // 35
  AliDielectronSignalMC* piKaonPairs = new AliDielectronSignalMC("piKaonPairs","pion+kaon pairs");  // mis-id pion + mis-id kaon
  piKaonPairs->SetLegPDGs(211,321);
  piKaonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(piKaonPairs);

  // 36
  AliDielectronSignalMC* piProtonPairs = new AliDielectronSignalMC("piProtonPairs","pion+proton pairs");  // mis-id pion + mis-id proton
  piProtonPairs->SetLegPDGs(211,2212);
  piProtonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(piProtonPairs);

  // 37
  AliDielectronSignalMC* kaonKaonPairs = new AliDielectronSignalMC("kaonKaonPairs","kaon+kaon pairs");  // mis-id kaon + mis-id kaon
  kaonKaonPairs->SetLegPDGs(321,321);
  kaonKaonPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(kaonKaonPairs);

  // 38
  AliDielectronSignalMC* muonAllPairs = new AliDielectronSignalMC("muonAllPairs","muon+everything pairs");        // mis-id muon + something else (electron, pion, kaon, proton)
  muonAllPairs->SetLegPDGs(13,13,kFALSE,kTRUE);
  muonAllPairs->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(muonAllPairs);
}

//______________________________________________________________________________________
void SetEtaCorrection()
{
//
// Eta equalization for the TPC response
//
  if (AliDielectronPID::GetEtaCorrFunction()) return;

  TString list=gSystem->Getenv("LIST");

  TFile f("$TRAIN_ROOT/jpsi_JPSI/EtaCorrMaps.root");
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
