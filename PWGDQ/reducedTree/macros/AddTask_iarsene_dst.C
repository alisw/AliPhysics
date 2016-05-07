void AddFMDTask();

//__________________________________________________________________________________________
AliAnalysisTask *AddTask_iarsene_dst(Int_t reducedEventType=-1, Bool_t writeTree=kTRUE){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_iarsene_dst", "No analysis manager found.");
    return 0;
  }

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  
  //create task and add it to the manager
  AliAnalysisTaskReducedTreeMaker *task=new AliAnalysisTaskReducedTreeMaker("DSTTreeMaker", kTRUE);
  //task->SetTriggerMask(AliVEvent::kMB);
  task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral);
  
  //if(trainConfig=="pp") task->SetRejectPileup();
  
  //task->SelectCollisionCandidates(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral);  // Events passing trigger and physics selection for analysis
  
  //task->UsePhysicsSelection(kTRUE);
  task->SetUseAnalysisUtils(kTRUE);
  
  //task->SetFillV0Info(kFALSE);
  //task->SetFillGammaConversions(kFALSE);
  //task->SetFillK0s(kFALSE);
  //task->SetFillLambda(kFALSE);
  //task->SetFillALambda(kFALSE);
  //task->SetFillCaloClusterInfo(kFALSE);
  //task->SetFillDielectronInfo(kFALSE);
  //task->SetFillFriendInfo(kFALSE);
  
  task->SetEventFilter(CreateEventFilter(isAOD));
  task->SetTrackFilter(CreateGlobalTrackFilter(isAOD));
  //  task->SetFlowTrackFilter(CreateFlowTrackFilter(isAOD));
  task->SetK0sPionCuts(CreateK0sPionCuts(isAOD));
  task->SetLambdaProtonCuts(CreateLambdaProtonCuts(isAOD));
  task->SetLambdaPionCuts(CreateLambdaPionCuts(isAOD));
  task->SetGammaElectronCuts(CreateGammaConvElectronCuts(isAOD));
  task->SetK0sMassRange(0.44,0.55);
  task->SetLambdaMassRange(1.090,1.14);
  task->SetGammaConvMassRange(0.0,0.1);
  task->SetV0OpenCuts(CreateV0OpenCuts(AliESDv0KineCuts::kPurity, AliESDv0KineCuts::kPbPb));
  task->SetV0StrongCuts(CreateV0StrongCuts(AliESDv0KineCuts::kPurity, AliESDv0KineCuts::kPbPb));
  //task->SetFillFMDInfo(); AddFMDTask();
  //task->SetFillMCInfo(kTRUE);
  //task->SetFillEventPlaneInfo(kTRUE);

  //task->SetTreeWritingOption(AliAnalysisTaskReducedTreeMaker::kFullEventsWithFullTracks);
  //task->SetTreeWritingOption(AliAnalysisTaskReducedTreeMaker::kBaseEventsWithBaseTracks);
  if(reducedEventType!=-1)
    task->SetTreeWritingOption(reducedEventType);
  task->SetWriteTree(writeTree);

  // You can define the active branches of the tree, all other branches will be set to inactive
  //task->SetTreeActiveBranch("Event");
  //task->SetTreeActiveBranch("fRunNo");
  // Alternatively, you can define the inactive branches of the tree, all other branches will be set to active
  //task->SetTreeInactiveBranch("fEventTag");
  //task->SetTreeInactiveBranch("fRunNo");
  //task->SetTreeInactiveBranch("fCandidates.*");
  
  mgr->AddTask(task);
    
  AliAnalysisDataContainer *cReducedTree = 0x0;
  if(task->WriteTree())
    cReducedTree = mgr->CreateContainer("dstTree", TTree::Class(),
                                          AliAnalysisManager::kOutputContainer, "dstTree.root");
  
  AliAnalysisDataContainer *cReducedEvent =
  mgr->CreateContainer("ReducedEventDQ",
                         AliReducedBaseEvent::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "reducedEvent");
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cReducedEvent);
  if(task->WriteTree()) mgr->ConnectOutput(task, 2, cReducedTree );
  
  return task;
}


//_______________________________________________________________________________________________
AliAnalysisCuts* CreateEventFilter(Bool_t isAOD) {
  //
  // Event wise cuts
  //
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  return eventCuts;
}


//_______________________________________________________________________________________________
AliAnalysisCuts* CreateGlobalTrackFilter(Bool_t isAOD) {
  //
  // Cuts for tracks to be written in the dst
  //
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  // general cuts ---------------------------------
  //if(isAOD) {
  AliDielectronVarCuts *trackCuts = new AliDielectronVarCuts("trackCuts","track cuts");
  trackCuts->AddCut(AliDielectronVarManager::kImpactParXY,-3.0,3.0);
  trackCuts->AddCut(AliDielectronVarManager::kImpactParZ,-10.0,10.0);
  trackCuts->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,0.15,1.0e+30);
  trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,1.0,161.0);
  //  trackCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.1,4.0);
  
  //trackCuts->AddCut(AliDielectronVarManager::kP,1.0,1.0e+30);
  //trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,70.0,161.0);
  cuts->AddCut(trackCuts);
    //}
  /*else {
    AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;
    // basic track quality cuts  (basicQ)
    esdTrackCuts->SetMaxDCAToVertexZ(10.0);
    esdTrackCuts->SetMaxDCAToVertexXY(10.0);
    esdTrackCuts->SetEtaRange(-0.90, 0.90);
    //  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    //esdTrackCuts->SetRequireITSRefit(kTRUE);
    //esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetPRange(0.2,1e30);
    //esdTrackCuts->SetPtRange(0.15,1e30);
    esdTrackCuts->SetMinNClustersTPC(50);
    esdTrackCuts->SetMaxChi2PerClusterTPC(4);
    //  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    cuts->AddCut(esdTrackCuts);
  }*/
  
  AliDielectronPID *electronPid = new AliDielectronPID("PID","PID cut");
  electronPid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-4.0, 4.0, 0.0, 0.0, kFALSE, AliDielectronPID::kRequire); // TPC 3-sigma inclusion for electron     
  electronPid->AddCut(AliDielectronPID::kTOF,AliPID::kProton,  -3.0, 3.0, -2.0, 2.0, kTRUE, AliDielectronPID::kRequire,AliDielectronVarManager::kTPCnSigmaPro); // TPC exclusion for proton   
  //cuts->AddCut(electronPid);
  
  return cuts;
}

//_________________________________________________________________________________________________________________
void AddFMDTask(){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  gSystem->Load("libPWGLFforward2");  // for FMD

  // Create the FMD task and add it to the manager
  //===========================================================================


  //--- AOD output handler -----------------------------------------
  AliAODHandler* ret = new AliAODHandler();
  //ret->SetFillAOD(kTRUE);
  ret->SetOutputFileName("AliAOD.pass2.root");
  mgr->SetOutputEventHandler(ret);
  ////AliVEventHandler*  outputHandlernew AliAODInputHandler();
  //


  //gROOT->LoadClass("AliAODForwardMult", "libPWGLFforward2");
  gSystem->Load("libESDfilter.so");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/ESDfilter/macros/AddTaskESDFilter.C");
  AliAnalysisTaskESDfilter *esdfilter = AddTaskESDFilter(kTRUE, kFALSE, kFALSE, kFALSE, kFALSE, kTRUE);
  //AliAnalysisTaskESDfilter *esdfilter = AddTaskESDFilter();

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (esdfilter,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (esdfilter,  0, mgr->GetCommonOutputContainer());

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/AddTaskForwardMult.C");

  Bool_t   mc  = false; // false: real data, true: simulated data
  ULong_t run = 0; // 0: get from data???
  UShort_t sys = 0; // 0: get from data, 1: pp, 2: AA 
  UShort_t sNN = 0; // 0: get from data, otherwise center of mass energy (per nucleon pair)
  Short_t  fld = 0; // 0: get from data, otherwise L3 field in kG
  //AliAnalysisTask *taskFmd  = AddTaskForwardMult(mc, sys, sNN, fld);
  //const Char_t* config = "$ALICE_PHYSICS/PWGDQ/dielectron/ReducedTreeMaker/ForwardAODConfig2.C";
  const Char_t* config = "ForwardAODConfig2.C";  // to be changed to the above line when config is committed
  AliAnalysisTask *taskFmd  = AddTaskForwardMult(mc, run, sys, sNN, fld, config);



  // --- Make the output container and connect it --------------------
  AliAnalysisDataContainer* histOut = 
    mgr->CreateContainer("Forward", TList::Class(), 
        AliAnalysisManager::kExchangeContainer,
        "Forward");
  //AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *output = 
    mgr->CreateContainer("ForwardResults", TList::Class(), 
        AliAnalysisManager::kParamContainer, 
        "ForwardResults");
  //AliAnalysisManager::GetCommonFileName());

  mgr->ConnectInput(taskFmd, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskFmd, 1, histOut);
  //mgr->ConnectOutput(taskFmd, 2, output);
  mgr->AddTask(taskFmd);
}

//_______________________________________________________________________________________________
AliAnalysisCuts* CreateFlowTrackFilter(Bool_t isAOD) {
  //
  // Cuts for tracks to be used for the event plane q-vector
  // These cuts are applied aditionally to the global track cuts
  //
  if(!isAOD) {
    AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts;
    esdTrackCuts->SetPtRange(0.2,2.0);
    esdTrackCuts->SetEtaRange(-0.8, 0.8);
    esdTrackCuts->SetMinNClustersTPC(70);
    return esdTrackCuts;
  }
  else return 0;
}


//_______________________________________________________________________________________________
AliAnalysisCuts* CreateK0sPionCuts(Bool_t isAOD) {
  //
  // Cuts on the K0s pions (tracking cuts, pid cuts, ...) 
  //
  if(!isAOD) {
    AliESDtrackCuts *pionCuts = new AliESDtrackCuts;
    pionCuts->SetPtRange(0.15,100.0);
    pionCuts->SetEtaRange(-1.6,1.6);
    pionCuts->SetRequireTPCRefit(kTRUE);
    pionCuts->SetMinNClustersTPC(50);
    return pionCuts;
  }
  else return 0;
}

//_______________________________________________________________________________________________
AliAnalysisCuts* CreateLambdaPionCuts(Bool_t isAOD) {
  //
  // Cuts on the Lambda pions (tracking cuts, pid cuts, ...) 
  //
  if(!isAOD) {
    AliESDtrackCuts *pionCuts = new AliESDtrackCuts;
    pionCuts->SetPtRange(0.15,100.0);
    pionCuts->SetEtaRange(-1.6,1.6);
    pionCuts->SetRequireTPCRefit(kTRUE);
    pionCuts->SetMinNClustersTPC(50);
    return pionCuts;
  }
  else return 0;
}


//_______________________________________________________________________________________________
AliAnalysisCuts* CreateLambdaProtonCuts(Bool_t isAOD) {
  //
  // Cuts on the Lambda protons (tracking cuts, pid cuts, ...) 
  //
  if(!isAOD) {
    AliESDtrackCuts *protonCuts = new AliESDtrackCuts;
    protonCuts->SetPtRange(0.15,100.0);
    protonCuts->SetEtaRange(-1.6,1.6);
    protonCuts->SetRequireTPCRefit(kTRUE);
    protonCuts->SetMinNClustersTPC(50);
    return protonCuts;
  }
  else return 0;
}


//______________________________________________________________________________________
AliAnalysisCuts* CreateGammaConvElectronCuts(Bool_t isAOD) {
  //
  // Cuts for the selection of electrons from gamma conversions
  //
  if(!isAOD) {
    AliESDtrackCuts* electronCuts = new AliESDtrackCuts;
    electronCuts->SetPtRange(0.15,100.0);
    electronCuts->SetEtaRange(-1.6,1.6);
    electronCuts->SetRequireTPCRefit(kTRUE);
    electronCuts->SetMinNClustersTPC(50);
    return electronCuts;
  }
}

//______________________________________________________________________________________
AliESDv0KineCuts* CreateV0StrongCuts(Int_t mode, Int_t type) {
  //
  // cuts for the selection of gamma conversions
  //
  AliESDv0KineCuts* cuts = new AliESDv0KineCuts();
  cuts->SetMode(mode, type);
  
  // leg cuts
  cuts->SetNTPCclusters(50);
  cuts->SetTPCrefit(kTRUE);
  cuts->SetTPCchi2perCls(4.0);
  cuts->SetTPCclusterratio(0.6);
  cuts->SetNoKinks(kTRUE);
  // gamma cuts                                                                                                                      
  cuts->SetGammaCutChi2NDF(10.0);
  Float_t cosPoint[2] = {0.0, 0.02};
  cuts->SetGammaCutCosPoint(cosPoint);
  Float_t cutDCA[2] = {0.0, 0.25};
  cuts->SetGammaCutDCA(cutDCA);
  Float_t vtxR[2] = {3.0, 90.0};
  cuts->SetGammaCutVertexR(vtxR);
  Float_t psiPairCut[2]={0.0,0.05};
  cuts->SetGammaCutPsiPair(psiPairCut);
  cuts->SetGammaCutInvMass(0.05);
  // K0s cuts
  cuts->SetK0CutChi2NDF(10.0);
  Float_t cosPointK0s[2] = {0.0, 0.02};
  cuts->SetK0CutCosPoint(cosPointK0s);
  Float_t cutDCAK0s[2] = {0.0, 0.2};
  cuts->SetK0CutDCA(cutDCAK0s);
  Float_t vtxRK0s[2] = {2.0, 30.0};
  cuts->SetK0CutVertexR(vtxRK0s);
  Float_t k0sInvMass[2] = {0.486, 0.508};
  cuts->SetK0CutInvMass(k0sInvMass);
  // Lambda and anti-Lambda cuts
  cuts->SetLambdaCutChi2NDF(10.0);
  Float_t cosPointLambda[2] = {0.0, 0.02};
  cuts->SetLambdaCutCosPoint(cosPointLambda);
  Float_t cutDCALambda[2] = {0.0, 0.2};
  cuts->SetLambdaCutDCA(cutDCALambda);
  Float_t vtxRLambda[2] = {2.0, 40.0};
  cuts->SetLambdaCutVertexR(vtxRLambda);
  Float_t lambdaInvMass[2] = {1.11, 1.12};
  cuts->SetLambdaCutInvMass(lambdaInvMass);
  
  return cuts;
}

//______________________________________________________________________________________
AliESDv0KineCuts* CreateV0OpenCuts(Int_t mode, Int_t type) {
  //
  // cuts for the selection of gamma conversions
  //
  AliESDv0KineCuts* cuts = new AliESDv0KineCuts();
  cuts->SetMode(mode, type);
  
  // leg cuts
  cuts->SetNTPCclusters(50);
  cuts->SetTPCrefit(kTRUE);
  cuts->SetTPCchi2perCls(4.0);
  cuts->SetTPCclusterratio(0.6);
  cuts->SetNoKinks(kTRUE);
  // gamma cuts                                                                                                                      
  cuts->SetGammaCutChi2NDF(10.0);
  Float_t cosPoint[2] = {0.0, 0.05};  //  Float_t cosPoint[2] = {0.0, 0.02};
  cuts->SetGammaCutCosPoint(cosPoint);
  Float_t cutDCA[2] = {0.0, 1.0};  //  Float_t cutDCA[2] = {0.0, 0.25};
  cuts->SetGammaCutDCA(cutDCA);
  Float_t vtxR[2] = {3.0, 90.0};
  cuts->SetGammaCutVertexR(vtxR);
  Float_t psiPairCut[2]={0.0,0.20};  //Float_t psiPairCut[2]={0.0,0.05};
  cuts->SetGammaCutPsiPair(psiPairCut);
  cuts->SetGammaCutInvMass(0.1); //  cuts->SetGammaCutInvMass(0.05);
  // K0s cuts
  cuts->SetK0CutChi2NDF(10.0);
  Float_t cosPointK0s[2] = {0.0, 0.1};  //  Float_t cosPointK0s[2] = {0.0, 0.02};
  cuts->SetK0CutCosPoint(cosPointK0s);
  Float_t cutDCAK0s[2] = {0.0, 1.0};  //  Float_t cutDCAK0s[2] = {0.0, 0.2};
  cuts->SetK0CutDCA(cutDCAK0s);
  Float_t vtxRK0s[2] = {2.0, 90.0};  //  Float_t vtxRK0s[2] = {2.0, 30.0};
  cuts->SetK0CutVertexR(vtxRK0s);
  Float_t k0sInvMass[2] = {0.44, 0.55};  //  Float_t k0sInvMass[2] = {0.486, 0.508};
  cuts->SetK0CutInvMass(k0sInvMass);
  // Lambda and anti-Lambda cuts
  cuts->SetLambdaCutChi2NDF(10.0);
  Float_t cosPointLambda[2] = {0.0, 0.1};  //  Float_t cosPointLambda[2] = {0.0, 0.02};
  cuts->SetLambdaCutCosPoint(cosPointLambda);
  Float_t cutDCALambda[2] = {0.0, 1.0};  //  Float_t cutDCALambda[2] = {0.0, 0.2};
  cuts->SetLambdaCutDCA(cutDCALambda);
  Float_t vtxRLambda[2] = {2.0, 90.0};  //  Float_t vtxRLambda[2] = {2.0, 40.0};
  cuts->SetLambdaCutVertexR(vtxRLambda);
  Float_t lambdaInvMass[2] = {1.09, 1.14};  //  Float_t lambdaInvMass[2] = {1.11, 1.12};
  cuts->SetLambdaCutInvMass(lambdaInvMass);
  
  return cuts;
}
