//===============================================================================================================
// addtask to create trees for J/psi - cross section analysis in pp 13TeV (last updated: 28/02/2019)
//===============================================================================================================

AliAnalysisCuts* EventFilter(Bool_t isAOD);
AliAnalysisCuts* JPsiElectronFilter(Bool_t isAOD);
AliAnalysisCuts* JPsiElectronFilterOpenPID(Bool_t isAOD);
AliAnalysisCuts* JPsiElectronFilterNoPID(Bool_t isAOD);
AliAnalysisCuts* JPsiElectronPreFilter(Bool_t isAOD);
AliAnalysisCuts* CreateK0sPionCuts(Bool_t isAOD);
AliAnalysisCuts* CreateLambdaPionCuts(Bool_t isAOD);
AliAnalysisCuts* CreateLambdaProtonCuts(Bool_t isAOD);
AliAnalysisCuts* CreateGammaConvElectronCuts(Bool_t isAOD);
AliESDv0KineCuts* CreateV0OpenCuts(Int_t mode, Int_t type);
AliESDv0KineCuts* CreateV0StrongCuts(Int_t mode, Int_t type);
void AddMCSignals(AliAnalysisTaskReducedTreeMaker* task);
void SetV0SelectionCuts(AliAnalysisTaskReducedTreeMaker *task, Bool_t isAOD);
void SetInactiveBranches(AliAnalysisTaskReducedTreeMaker *task);

//_______________________________________________________________________________________________________________
AliAnalysisTask *AddTask_ilofnes_dst(Int_t reducedEventType=-1, Bool_t writeTree=kTRUE, TString  prod="") {

    //get current analysis manager
    //-----------------------------------------------------------------------------------------------------------
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) { Error("AddTask_ilofnes_dst", "No analysis manager found."); return 0; }

    // query MC handler and AOD
    //-----------------------------------------------------------------------------------------------------------
    Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
    Bool_t isAOD = mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

    // select trigger according to production
    //-----------------------------------------------------------------------------------------------------------
    Int_t triggerChoice = 2; //Chose only kINT7
    //Int_t triggerChoice = 9;
    printf("AddTask_ilofnes_dst() trigger choice set to %d (%s)\n", triggerChoice, prod.Data());

    //create task and add it to the manager
    //-----------------------------------------------------------------------------------------------------------
    AliAnalysisTaskReducedTreeMaker *task=new AliAnalysisTaskReducedTreeMaker("DSTTreeMaker", kTRUE);
    if (triggerChoice==0)       task->SetTriggerMask(AliVEvent::kAny);
    else if (triggerChoice==1)  task->SetTriggerMask(AliVEvent::kMB);
    else if (triggerChoice==2)  task->SetTriggerMask(AliVEvent::kINT7);
    else if (triggerChoice==3)  task->SetTriggerMask(AliVEvent::kTRD);
    else if (triggerChoice==4)  task->SetTriggerMask(AliVEvent::kEMCEJE+AliVEvent::kEMCEGA);
    else if (triggerChoice==5)  task->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kTRD);
    else if (triggerChoice==6)  task->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kEMCEJE+AliVEvent::kEMCEGA);
    else if (triggerChoice==7)  task->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kTRD+AliVEvent::kEMCEJE+AliVEvent::kEMCEGA);
    else if (triggerChoice==8)  task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kINT7+AliVEvent::kHighMultSPD+AliVEvent::kHighMultV0);
    else if (triggerChoice==9)  task->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kTRD+AliVEvent::kEMCEJE+AliVEvent::kEMCEGA+AliVEvent::kHighMultSPD+AliVEvent::kHighMultV0);
    else { printf("WARNING: In AddTask_ilofnes_dst(), no trigger specified, or not supported! Using kINT7.\n"); task->SetTriggerMask(AliVEvent::kINT7); }

    // pile-up, physics selection and analysis utils
    //-----------------------------------------------------------------------------------------------------------
    //task->SetRejectPileup(kTRUE);
    task->UsePhysicsSelection(kTRUE); //not if MC??
    task->SetUseAnalysisUtils(kTRUE);

    // toggle filling of branches of tree
    //-----------------------------------------------------------------------------------------------------------
    task->SetFillTrackInfo(kTRUE);
    task->SetFillV0Info(kTRUE);  //Set to kTRUE??
    task->SetFillGammaConversions(kTRUE);
    task->SetFillK0s(kTRUE);
    task->SetFillLambda(kTRUE);
    task->SetFillALambda(kTRUE);
    task->SetFillCaloClusterInfo(kFALSE);
    task->SetFillFMDInfo(kFALSE);
    task->SetFillEventPlaneInfo(kFALSE,1.0); //tpcEtaGap = 1.0 (default)
    task->SetFillTRDMatchedTracks(kFALSE);
    if (hasMC) {
      task->SetFillMCInfo(kTRUE);
      AddMCSignals(task);
    }

    // event selection
    //-----------------------------------------------------------------------------------------------------------
    task->SetEventFilter(EventFilter(isAOD));

    // track selection
    //-----------------------------------------------------------------------------------------------------------
    if (hasMC)
      task->AddTrackFilter(JPsiElectronFilterNoPID(isAOD),  kFALSE,2);  // full track info, fQualityFlag bit 32 - tracking cut, no PID cuts
    else  {
      task->AddTrackFilter(JPsiElectronFilter(isAOD),       kFALSE,2);  // full track info, fQualityFlag bit 32
      task->AddTrackFilter(JPsiElectronFilterOpenPID(isAOD),  kFALSE,2);  // full track info, open PID, fQualityFlag bit 33
    }
    task->AddTrackFilter(JPsiElectronPreFilter(isAOD),    kTRUE);  // full track info, fQualityFlag bit 34 - prefilter

    // V0 selection
    //-----------------------------------------------------------------------------------------------------------
    SetV0SelectionCuts(task, isAOD);

    // active/inactive branches of tree
    //-----------------------------------------------------------------------------------------------------------
    SetInactiveBranches(task);

    // set writing options
    //  AliAnalysisTaskReducedTreeMaker::kBaseEventsWithBaseTracks = 0
    //  AliAnalysisTaskReducedTreeMaker::kBaseEventsWithFullTracks = 1
    //  AliAnalysisTaskReducedTreeMaker::kFullEventsWithBaseTracks = 2
    //  AliAnalysisTaskReducedTreeMaker::kFullEventsWithFullTracks = 3
    //-----------------------------------------------------------------------------------------------------------
    if(reducedEventType!=-1) task->SetTreeWritingOption(reducedEventType);
    task->SetWriteTree(writeTree);
    task->SetWriteUnbiasedEvents(0.002);

    // add task to manager
    //-----------------------------------------------------------------------------------------------------------
    mgr->AddTask(task);

    // connect output containers
    //-----------------------------------------------------------------------------------------------------------
    AliAnalysisDataContainer *cReducedTree    = 0x0;
    AliAnalysisDataContainer *cEventStatsInfo = 0x0;

    if(task->WriteTree()) {
      cReducedTree            = mgr->CreateContainer("dstTree",       TTree::Class(),AliAnalysisManager::kOutputContainer, "dstTree.root");
      cEventStatsInfo         = mgr->CreateContainer("EventStats",    TList::Class(), AliAnalysisManager::kOutputContainer, "dstTree.root");
    }
    AliAnalysisDataContainer *cReducedEvent = mgr->CreateContainer("ReducedEventDQ", AliReducedBaseEvent::Class(), AliAnalysisManager::kExchangeContainer, "reducedEvent");
    mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, cReducedEvent);
    if(task->WriteTree()) {
      mgr->ConnectOutput(task, 2, cReducedTree);
      mgr->ConnectOutput(task, 3, cEventStatsInfo);
    }

    // done
    //-----------------------------------------------------------------------------------------------------------
    return task;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* EventFilter(Bool_t isAOD) {
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  //eventCuts->SetVertexZ(-10.,10.);
  return eventCuts;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* JPsiElectronFilter(Bool_t isAOD) {
  AliDielectronCutGroup*  jpsiElectrons = new AliDielectronCutGroup("jpsiElectrons","J/psi candidate electrons (basic cut)",AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts*   trackCuts     = new AliDielectronVarCuts("trackCuts0","track cuts");
  trackCuts->AddCut(AliDielectronVarManager::kImpactParXY,-1.0,1.0);
  trackCuts->AddCut(AliDielectronVarManager::kImpactParZ,-3.0,3.0);
  trackCuts->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,1.0,1.0e+30);
  trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,70.0,161.0);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaEle, -4.0, 4.0);
  trackCuts->AddCut(AliDielectronVarManager::kPIn,0.0,5.0);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPro, -2000.0, 2.0, kTRUE);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPio, -2000.0, 2.0, kTRUE);
  jpsiElectrons->AddCut(trackCuts);
  AliDielectronTrackCuts* trackCuts2 = new AliDielectronTrackCuts("trackCuts1","track cuts");
  trackCuts2->SetRequireITSRefit(kTRUE);
  trackCuts2->SetRequireTPCRefit(kTRUE);
  if(isAOD) trackCuts2->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  jpsiElectrons->AddCut(trackCuts2);
  return jpsiElectrons;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* JPsiElectronPreFilter(Bool_t isAOD) {
  AliDielectronCutGroup* jpsiElectrons = new AliDielectronCutGroup("jpsiElectronsPrefilterCut","J/psi candidate electrons(prefilter cut)",AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts *trackCuts = new AliDielectronVarCuts("trackCuts2","track cuts");
  trackCuts->AddCut(AliDielectronVarManager::kImpactParXY,-3.0,3.0);
  trackCuts->AddCut(AliDielectronVarManager::kImpactParZ,-10.0,10.0);
  trackCuts->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,0.3,1.0e+30);
  trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,50.0,161.0);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaEle, -4.0, 4.0);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPro, -2000.0, 2.0, kTRUE);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPio, -2000.0, 2.0, kTRUE);
  jpsiElectrons->AddCut(trackCuts);
  AliDielectronTrackCuts* trackCuts2 = new AliDielectronTrackCuts("trackCuts3","track cuts");
  //trackCuts2->SetRequireITSRefit(kTRUE);
  trackCuts2->SetRequireTPCRefit(kTRUE);
  if(isAOD) trackCuts2->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  jpsiElectrons->AddCut(trackCuts2);
  return jpsiElectrons;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* JPsiElectronFilterOpenPID(Bool_t isAOD) {
  AliDielectronCutGroup*  jpsiElectrons = new AliDielectronCutGroup("jpsiElectronsOpenPID","J/psi candidate electrons (basic cut)",AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts*   trackCuts     = new AliDielectronVarCuts("trackCuts0","track cuts 0");
  trackCuts->AddCut(AliDielectronVarManager::kImpactParXY,-1.0,1.0);
  trackCuts->AddCut(AliDielectronVarManager::kImpactParZ,-3.0,3.0);
  trackCuts->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,1.0,1.0e+30);
  trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,70.0,161.0);
  trackCuts->AddCut(AliDielectronVarManager::kPIn,5.0,1.0e+30);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaEle, -4.0, 4.0);
  jpsiElectrons->AddCut(trackCuts);
  AliDielectronTrackCuts* trackCuts2 = new AliDielectronTrackCuts("trackCuts1","track cuts 1");
  trackCuts2->SetRequireITSRefit(kTRUE);
  trackCuts2->SetRequireTPCRefit(kTRUE);
  if(isAOD) trackCuts2->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  jpsiElectrons->AddCut(trackCuts2);
  return jpsiElectrons;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* JPsiElectronFilterNoPID(Bool_t isAOD) {
  AliDielectronCutGroup* jpsiElectrons = new AliDielectronCutGroup("jpsiElectronsNoPID","J/psi candidate electrons (No PID)",AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts *trackCuts = new AliDielectronVarCuts("trackCuts4","track cuts");
  trackCuts->AddCut(AliDielectronVarManager::kImpactParXY,-1.0,1.0);
  trackCuts->AddCut(AliDielectronVarManager::kImpactParZ,-3.0,3.0);
  trackCuts->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,1.0,1.0e+30);
  trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,70.0,161.0);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaEle, -5.0,5.0);
  jpsiElectrons->AddCut(trackCuts);
  AliDielectronTrackCuts* trackCuts2 = new AliDielectronTrackCuts("trackCuts5","track cuts");
  trackCuts2->SetRequireITSRefit(kTRUE);
  trackCuts2->SetRequireTPCRefit(kTRUE);
  if(isAOD) trackCuts2->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  jpsiElectrons->AddCut(trackCuts2);
  return jpsiElectrons;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* CreateK0sPionCuts(Bool_t isAOD) {
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

//_______________________________________________________________________________________________________________
AliAnalysisCuts* CreateLambdaPionCuts(Bool_t isAOD) {
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

//_______________________________________________________________________________________________________________
AliAnalysisCuts* CreateLambdaProtonCuts(Bool_t isAOD) {
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

//_______________________________________________________________________________________________________________
AliAnalysisCuts* CreateGammaConvElectronCuts(Bool_t isAOD) {
    AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cutsEleV0","cuts",AliDielectronCutGroup::kCompAND);
    if(!isAOD) {
        AliESDtrackCuts* electronCuts = new AliESDtrackCuts;
        electronCuts->SetPtRange(0.15,100.0);
        electronCuts->SetEtaRange(-1.6,1.6);
        electronCuts->SetRequireTPCRefit(kTRUE);
        electronCuts->SetMinNClustersTPC(50);
        cuts->AddCut(electronCuts);
    }
    AliDielectronPID *electronPid = new AliDielectronPID("PIDeleV0","PID cut");
    electronPid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-4.0, 4.0, 0.0, 0.0, kFALSE, AliDielectronPID::kRequire); // TPC 4-sigma inclusion for electron
    cuts->AddCut(electronPid);
    return cuts;
}

//_______________________________________________________________________________________________________________
AliESDv0KineCuts* CreateV0OpenCuts(Int_t mode, Int_t type) {
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
    Float_t cosPoint[2] = {0.0, 0.05};
    cuts->SetGammaCutCosPoint(cosPoint);
    Float_t cutDCA[2] = {0.0, 1.0};
    cuts->SetGammaCutDCA(cutDCA);
    Float_t vtxR[2] = {3.0, 90.0};
    cuts->SetGammaCutVertexR(vtxR);
    Float_t psiPairCut[2]={0.0,0.20};
    cuts->SetGammaCutPsiPair(psiPairCut);
    cuts->SetGammaCutInvMass(0.1);
    /*
    // K0s cuts
    cuts->SetK0CutChi2NDF(10.0);
    Float_t cosPointK0s[2] = {0.0, 0.1};
    cuts->SetK0CutCosPoint(cosPointK0s);
    Float_t cutDCAK0s[2] = {0.0, 1.0};
    cuts->SetK0CutDCA(cutDCAK0s);
    Float_t vtxRK0s[2] = {2.0, 90.0};
    cuts->SetK0CutVertexR(vtxRK0s);
    Float_t k0sInvMass[2] = {0.44, 0.55};
    cuts->SetK0CutInvMass(k0sInvMass);
    // Lambda and anti-Lambda cuts
    cuts->SetLambdaCutChi2NDF(10.0);
    Float_t cosPointLambda[2] = {0.0, 0.1};
    cuts->SetLambdaCutCosPoint(cosPointLambda);
    Float_t cutDCALambda[2] = {0.0, 1.0};
    cuts->SetLambdaCutDCA(cutDCALambda);
    Float_t vtxRLambda[2] = {2.0, 90.0};
    cuts->SetLambdaCutVertexR(vtxRLambda);
    Float_t lambdaInvMass[2] = {1.09, 1.14};
    cuts->SetLambdaCutInvMass(lambdaInvMass);
     */
    return cuts;
}

//_______________________________________________________________________________________________________________
AliESDv0KineCuts* CreateV0StrongCuts(Int_t mode, Int_t type) {
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
      /*
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
       */
    return cuts;
}

//_______________________________________________________________________________________________________________
void AddMCSignals(AliAnalysisTaskReducedTreeMaker* task) {
  //
  // Add the MC signals to be filtered
  //

  // 0 = kJpsiInclusive
  AliSignalMC* jpsiInclusive = new AliSignalMC("JpsiInclusive", "",1,1);
  jpsiInclusive->SetPDGcode(0, 0, 443, kFALSE);
  task->AddMCsignal(jpsiInclusive, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 1 = kJpsiNonPrompt
  AliSignalMC* jpsiFromB = new AliSignalMC("JpsiNonPrompt","",1,2);
  jpsiFromB->SetPDGcode(0, 0, 443, kFALSE);
  jpsiFromB->SetPDGcode(0, 1, 503, kTRUE); //503 - all beauty hadrons
  task->AddMCsignal(jpsiFromB, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 2 = kJpsiPrompt
  AliSignalMC* jpsiPrompt = new AliSignalMC("JpsiPrompt","",1,2);
  jpsiPrompt->SetPDGcode(0, 0, 443, kFALSE);
  jpsiPrompt->SetPDGcode(0, 1, 503, kTRUE, kTRUE);
  task->AddMCsignal(jpsiPrompt, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 3 = kJpsiRadiative
  AliSignalMC* jpsiInclusiveRadiative = new AliSignalMC("JpsiInclusiveRadiative","",1,1);
  jpsiInclusiveRadiative->SetPDGcode(0, 0, 443, kFALSE);
  jpsiInclusiveRadiative->SetSourceBit(0, 0, AliSignalMC::kRadiativeDecay);
  task->AddMCsignal(jpsiInclusiveRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 4 = kJpsiNonRadiative
  AliSignalMC* jpsiInclusiveNonRadiative = new AliSignalMC("JpsiInclusiveNonRadiative","",1,1);
  jpsiInclusiveNonRadiative->SetPDGcode(0, 0, 443, kFALSE);
  jpsiInclusiveNonRadiative->SetSourceBit(0, 0, AliSignalMC::kRadiativeDecay, kTRUE);
  task->AddMCsignal(jpsiInclusiveNonRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 5 = kJpsiNonPromptRadiative
  AliSignalMC* jpsiFromBRadiative = new AliSignalMC("JpsiFromBRadiative","",1,2);
  jpsiFromBRadiative->SetPDGcode(0, 0, 443, kFALSE);
  jpsiFromBRadiative->SetPDGcode(0, 1, 503, kTRUE);
  jpsiFromBRadiative->SetSourceBit(0, 0, AliSignalMC::kRadiativeDecay, kFALSE);
  task->AddMCsignal(jpsiFromBRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 6 = kJpsiNonPromptNonRadiative
  AliSignalMC* jpsiFromBNonRadiative = new AliSignalMC("JpsiFromBNonRadiative","",1,2);
  jpsiFromBNonRadiative->SetPDGcode(0, 0, 443, kFALSE);
  jpsiFromBNonRadiative->SetPDGcode(0, 1, 503, kTRUE);
  jpsiFromBNonRadiative->SetSourceBit(0, 0, AliSignalMC::kRadiativeDecay, kTRUE);
  task->AddMCsignal(jpsiFromBNonRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 7 = kJpsiDecayElectron
  AliSignalMC* electronFromJpsi = new AliSignalMC("electronFromJpsiInclusive","",1,2);
  electronFromJpsi->SetPDGcode(0, 0, 11, kTRUE);
  electronFromJpsi->SetPDGcode(0, 1, 443);
  task->AddMCsignal(electronFromJpsi, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 8 = kJpsiNonPromptDecayElectron
  AliSignalMC* electronFromJpsiNonPrompt = new AliSignalMC("electronFromJpsiNonPrompt","",1,3);
  electronFromJpsiNonPrompt->SetPDGcode(0, 0, 11, kTRUE);
  electronFromJpsiNonPrompt->SetPDGcode(0, 1, 443);
  electronFromJpsiNonPrompt->SetPDGcode(0, 2, 503, kTRUE);
  task->AddMCsignal(electronFromJpsiNonPrompt, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 9 = kJpsiPromptDecayElectron
  AliSignalMC* electronFromJpsiPrompt = new AliSignalMC("electronFromJpsiPrompt","", 1,3);
  electronFromJpsiPrompt->SetPDGcode(0, 0, 11, kTRUE);
  electronFromJpsiPrompt->SetPDGcode(0, 1, 443);
  electronFromJpsiPrompt->SetPDGcode(0, 2, 503, kTRUE, kTRUE);
  task->AddMCsignal(electronFromJpsiPrompt, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 10 = kJpsiRadiativeDecayElectron
  AliSignalMC* electronFromJpsiRadiative = new AliSignalMC("electronFromJpsiRadiative","",1,2);
  electronFromJpsiRadiative->SetPDGcode(0, 0, 11, kTRUE);
  electronFromJpsiRadiative->SetPDGcode(0, 1, 443);
  electronFromJpsiRadiative->SetSourceBit(0, 1, AliSignalMC::kRadiativeDecay);
  task->AddMCsignal(electronFromJpsiRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);
  // 11 = kJpsiNonRadiativeDecayElectron
  AliSignalMC* electronFromJpsiNonRadiative = new AliSignalMC("electronFromJpsiNonRadiative","",1,2);
  electronFromJpsiNonRadiative->SetPDGcode(0, 0, 11, kTRUE);
  electronFromJpsiNonRadiative->SetPDGcode(0, 1, 443);
  electronFromJpsiNonRadiative->SetSourceBit(0, 1, AliSignalMC::kRadiativeDecay, kTRUE);
  task->AddMCsignal(electronFromJpsiNonRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 12 = kJpsiDecayPhoton
  AliSignalMC* photonFromJpsiDecay = new AliSignalMC("photonFromJpsiDecay","",1,2);
  photonFromJpsiDecay->SetPDGcode(0, 0, 22);
  photonFromJpsiDecay->SetPDGcode(0, 1, 443);
  task->AddMCsignal(photonFromJpsiDecay, AliAnalysisTaskReducedTreeMaker::kFullTrack);
}

//_______________________________________________________________________________________________________________
void SetV0SelectionCuts(AliAnalysisTaskReducedTreeMaker *task, Bool_t isAOD) {
  task->SetK0sPionCuts(CreateK0sPionCuts(isAOD));
  task->SetLambdaPionCuts(CreateLambdaPionCuts(isAOD));
  task->SetLambdaProtonCuts(CreateLambdaProtonCuts(isAOD));
  task->SetGammaElectronCuts(CreateGammaConvElectronCuts(isAOD));
  task->SetK0sMassRange(0.44,0.55);
  task->SetLambdaMassRange(1.090,1.14);
  task->SetGammaConvMassRange(0.0,0.1);
  task->SetV0OpenCuts(CreateV0OpenCuts(AliESDv0KineCuts::kPurity, AliESDv0KineCuts::kPP));
  task->SetV0StrongCuts(CreateV0StrongCuts(AliESDv0KineCuts::kPurity, AliESDv0KineCuts::kPP));
}

//_______________________________________________________________________________________________________________
void SetInactiveBranches(AliAnalysisTaskReducedTreeMaker *task) {

  // event
  task->SetTreeInactiveBranch("fCentrality*");
  task->SetTreeInactiveBranch("fCentQuality");
  //task->SetTreeInactiveBranch("fNV0candidates*");
  //task->SetTreeInactiveBranch("fCandidates.*");
  task->SetTreeInactiveBranch("fEventNumberInFile");
  task->SetTreeInactiveBranch("fL0TriggerInputs");
  task->SetTreeInactiveBranch("fL1TriggerInputs");
  task->SetTreeInactiveBranch("fL2TriggerInputs");
  //task->SetTreeInactiveBranch("fBC");
  task->SetTreeInactiveBranch("fTimeStamp");
  task->SetTreeInactiveBranch("fEventType");
  task->SetTreeInactiveBranch("fMultiplicityEstimators*");
  task->SetTreeInactiveBranch("fMultiplicityEstimatorPercentiles*");
  task->SetTreeInactiveBranch("fIRIntClosestIntMap*");
  task->SetTreeInactiveBranch("fVtxTPC*");
  task->SetTreeInactiveBranch("fNVtxTPCContributors");
  task->SetTreeInactiveBranch("fVtxSPD*");
  task->SetTreeInactiveBranch("fNVtxSPDContributors");
  //task->SetTreeInactiveBranch("fNpileupSPD");
  //task->SetTreeInactiveBranch("fNpileupTracks");
  //task->SetTreeInactiveBranch("fNTPCclusters");
  task->SetTreeInactiveBranch("fNPMDtracks");
  task->SetTreeInactiveBranch("fNTRDtracks");
  task->SetTreeInactiveBranch("fNTRDtracklets");
  task->SetTreeInactiveBranch("fSPDntrackletsEta*");
  //task->SetTreeInactiveBranch("fITSClusters*");
  //task->SetTreeInactiveBranch("fNtracksPerTrackingFlag*");
  task->SetTreeInactiveBranch("fVZEROMult*");
  task->SetTreeInactiveBranch("fZDCnEnergy*");
  task->SetTreeInactiveBranch("fZDCpEnergy*");
  task->SetTreeInactiveBranch("fZDCnTotalEnergy*");
  task->SetTreeInactiveBranch("fZDCpTotalEnergy*");
  task->SetTreeInactiveBranch("fT0amplitude*");
  //task->SetTreeInactiveBranch("fT0TOF*");
  //task->SetTreeInactiveBranch("fT0TOFbest*");
  task->SetTreeInactiveBranch("fT0zVertex");
  //task->SetTreeInactiveBranch("fT0start");
  //task->SetTreeInactiveBranch("fT0pileup");
  task->SetTreeInactiveBranch("fT0sattelite");
  //task->SetTreeInactiveBranch("fNCaloClusters");
  //task->SetTreeInactiveBranch("fCaloClusters.*");
  task->SetTreeInactiveBranch("fFMD.*");
  task->SetTreeInactiveBranch("fEventPlane.*");

  // tracks
  task->SetTreeInactiveBranch("fTracks.fTPCPhi");
  task->SetTreeInactiveBranch("fTracks.fTPCPt");
  task->SetTreeInactiveBranch("fTracks.fTPCEta");
  //task->SetTreeInactiveBranch("fTracks.fMomentumInner");
  //task->SetTreeInactiveBranch("fTracks.fDCA*");
  task->SetTreeInactiveBranch("fTracks.fTPCDCA*");
  task->SetTreeInactiveBranch("fTracks.fTrackLength");
  task->SetTreeInactiveBranch("fTracks.fMassForTracking");
  //task->SetTreeInactiveBranch("fTracks.fChi2TPCConstrainedVsGlobal");
  task->SetTreeInactiveBranch("fTracks.fHelixCenter*");
  task->SetTreeInactiveBranch("fTracks.fHelixRadius");
  //task->SetTreeInactiveBranch("fTracks.fITSclusterMap");
  //task->SetTreeInactiveBranch("fTracks.fITSSharedClusterMap");
  task->SetTreeInactiveBranch("fTracks.fITSsignal");
  task->SetTreeInactiveBranch("fTracks.fITSnSig*");
  task->SetTreeInactiveBranch("fTracks.fTPCActiveLength");
  task->SetTreeInactiveBranch("fTracks.fTPCGeomLength");
  //task->SetTreeInactiveBranch("fTracks.fTOFbeta");
  //task->SetTreeInactiveBranch("fTracks.fTOFtime");
  //task->SetTreeInactiveBranch("fTracks.fTOFdx");
  //task->SetTreeInactiveBranch("fTracks.fTOFdz");
  //task->SetTreeInactiveBranch("fTracks.fTOFmismatchProbab");
  //task->SetTreeInactiveBranch("fTracks.fTOFchi2");
  //task->SetTreeInactiveBranch("fTracks.fTOFnSig*");
  //task->SetTreeInactiveBranch("fTracks.fTOFdeltaBC");
  task->SetTreeInactiveBranch("fTracks.fTRDntracklets*");
  task->SetTreeInactiveBranch("fTracks.fTRDpid*");
  task->SetTreeInactiveBranch("fTracks.fTRDpidLQ2D*");
  //task->SetTreeInactiveBranch("fTracks.fCaloClusterId");
  //task->SetTreeInactiveBranch("fTracks.fMCMom*");
  //task->SetTreeInactiveBranch("fTracks.fMCFreezeout*");
  //task->SetTreeInactiveBranch("fTracks.fMCLabels*");
  //task->SetTreeInactiveBranch("fTracks.fMCPdg*");
  //task->SetTreeInactiveBranch("fTracks.fMCGeneratorIndex");

  // candidates
  //task->SetTreeInactiveBranch("fCandidates.fP*");
  //task->SetTreeInactiveBranch("fCandidates.fIsCartesian");
  //task->SetTreeInactiveBranch("fCandidates.fCharge");
  //task->SetTreeInactiveBranch("fCandidates.fFlags");
  //task->SetTreeInactiveBranch("fCandidates.fQualityFlags");
  //task->SetTreeInactiveBranch("fCandidates.fCandidateId");
  //task->SetTreeInactiveBranch("fCandidates.fPairType");
  //task->SetTreeInactiveBranch("fCandidates.fLegIds*");
  //task->SetTreeInactiveBranch("fCandidates.fMass*");
  //task->SetTreeInactiveBranch("fCandidates.fLxy");
  //task->SetTreeInactiveBranch("fCandidates.fPointingAngle");
  //task->SetTreeInactiveBranch("fCandidates.fChisquare");
}
