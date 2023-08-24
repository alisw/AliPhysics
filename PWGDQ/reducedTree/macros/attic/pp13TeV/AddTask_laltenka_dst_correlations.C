//===============================================================================================================
// addtask to create trees for J/psi - hadron correlation analysis in pp 13TeV (last updated: 26/04/2019)
//===============================================================================================================

AliAnalysisCuts* EventFilter(Bool_t isAOD);
AliAnalysisCuts* CaloClusterFilter(Bool_t isAOD);
AliAnalysisCuts* JPsiElectronFilter(Bool_t isAOD);
AliAnalysisCuts* JPsiElectronPreFilter(Bool_t isAOD);
AliAnalysisCuts* AssociatedHadronFilter(Bool_t isAOD, Int_t spec);
AliAnalysisCuts* AssociatedPionFilter(Bool_t isAOD, Int_t spec);
void AddMCSignals(AliAnalysisTaskReducedTreeMaker* task);
void SetInactiveBranches(AliAnalysisTaskReducedTreeMaker *task);

//_______________________________________________________________________________________________________________
AliAnalysisTask *AddTask_laltenka_dst_correlations(Int_t reducedEventType=-1, Bool_t writeTree=kTRUE, TString  prod="") {

    //get current analysis manager
    //-----------------------------------------------------------------------------------------------------------
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) { Error("AddTask_laltenka_dst", "No analysis manager found."); return 0; }

    // query MC handler and AOD
    //-----------------------------------------------------------------------------------------------------------
    Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
    Bool_t isAOD = mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

    //create task
    //-----------------------------------------------------------------------------------------------------------
    AliAnalysisTaskReducedTreeMaker* task = new AliAnalysisTaskReducedTreeMaker("DSTTreeMaker", kTRUE);
  
    // select trigger (according to production)
    //-----------------------------------------------------------------------------------------------------------
    Int_t triggerChoice = 10;
    printf("AddTask_laltenka_dst() trigger choice set to %d (%s)\n", triggerChoice, prod.Data());
    if (triggerChoice==0)       // MB (INT7)
      task->SetTriggerMask(AliVEvent::kINT7);
    else if (triggerChoice==1)  // high mult. (SPD, V0)
      task->SetTriggerMask(AliVEvent::kHighMultSPD+AliVEvent::kHighMultV0);
    else if (triggerChoice==2)  // EMCal (L0, L1)
      task->SetTriggerMask(AliVEvent::kEMC7+AliVEvent::kEMCEGA);
    else if (triggerChoice==3)  // TRD
      task->SetTriggerMask(AliVEvent::kTRD);
    else if (triggerChoice==4)  // MB + high mult.
      task->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kHighMultSPD+AliVEvent::kHighMultV0);
    else if (triggerChoice==5)  // MB + EMCal
      task->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kEMC7+AliVEvent::kEMCEGA);
    else if (triggerChoice==6)  // MB + TRD
      task->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kTRD);
    else if (triggerChoice==9)  // MB + high mult. + EMCal + TRD
      task->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kHighMultSPD+AliVEvent::kHighMultV0+AliVEvent::kEMC7+AliVEvent::kEMCEGA+AliVEvent::kTRD);
    else if (triggerChoice==10)  // MB + high mult. + EMCal
      task->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kHighMultSPD+AliVEvent::kHighMultV0+AliVEvent::kEMCEGA);
    else {
      printf("WARNING: In AddTask_laltenka_dst(), no trigger specified, or not supported! Using kINT7.\n");
      task->SetTriggerMask(AliVEvent::kINT7);
    }

    // pile-up, physics selection and analysis utils
    //-----------------------------------------------------------------------------------------------------------
    task->SetRejectPileup(kFALSE);    // done at analysis level
    task->UsePhysicsSelection(kTRUE);
    task->SetUseAnalysisUtils(kTRUE);

    // toggle filling of branches of tree
    //-----------------------------------------------------------------------------------------------------------
    task->SetFillTrackInfo(kTRUE);
    task->SetFillCaloClusterInfo(kTRUE);
    task->SetFillV0Info(kFALSE);
    task->SetFillGammaConversions(kFALSE);
    task->SetFillK0s(kFALSE);
    task->SetFillLambda(kFALSE);
    task->SetFillALambda(kFALSE);
    task->SetFillFMDInfo(kFALSE);
    task->SetFillEventPlaneInfo(kFALSE);
    task->SetFillTRDMatchedTracks(kFALSE);
    if (hasMC) {
      task->SetFillMCInfo(kTRUE);
      AddMCSignals(task);
    }

    // event selection
    //-----------------------------------------------------------------------------------------------------------
    task->SetEventFilter(EventFilter(isAOD));

    // cluster selection
    //-----------------------------------------------------------------------------------------------------------
    task->AddCaloClusterFilter(CaloClusterFilter(isAOD));
  
    // track selection
    //-----------------------------------------------------------------------------------------------------------
    // electron tracks
    task->AddTrackFilter(JPsiElectronFilter(isAOD),         kFALSE, 2); // full track info, fQualityFlag bit 32 - low momentum
    task->AddTrackFilter(JPsiElectronPreFilter(isAOD),      kTRUE);     // base track info, fQualityFlag bit 33 - prefilter
    // associated unidentified hadrons
    task->AddTrackFilter(AssociatedHadronFilter(isAOD, 0),  kTRUE);     // base track info, fQualityFlag bit 34 - unid. hadron, ITS + TPC refit
    task->AddTrackFilter(AssociatedHadronFilter(isAOD, 1),  kTRUE);     // base track info, fQualityFlag bit 35 - unid. hadron, ITS + TPC refit, SPD any
    task->AddTrackFilter(AssociatedHadronFilter(isAOD, 2),  kTRUE);     // base track info, fQualityFlag bit 36 - unid. hadron, ITS + TPC refit, no kinks
    task->AddTrackFilter(AssociatedHadronFilter(isAOD, 3),  kTRUE);     // base track info, fQualityFlag bit 37 - unid. hadron, ITS + TPC refit, SPD any, no kinks
    // associated pions
    task->AddTrackFilter(AssociatedPionFilter(isAOD, 0),    kTRUE);     // base track info, fQualityFlag bit 38 - pions, ITS + TPC refit
    task->AddTrackFilter(AssociatedPionFilter(isAOD, 1),    kTRUE);     // base track info, fQualityFlag bit 39 - pions, ITS + TPC refit, SPD any
    task->AddTrackFilter(AssociatedPionFilter(isAOD, 2),    kTRUE);     // base track info, fQualityFlag bit 40 - pions, ITS + TPC refit, no kinks
    task->AddTrackFilter(AssociatedPionFilter(isAOD, 3),    kTRUE);     // base track info, fQualityFlag bit 41 - pions, ITS + TPC refit, SPD any, no kinks

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
    task->SetWriteUnbiasedEvents(0.02);

    // add task to manager
    //-----------------------------------------------------------------------------------------------------------
    mgr->AddTask(task);

    // connect output containers
    //-----------------------------------------------------------------------------------------------------------
    AliAnalysisDataContainer* cReducedEvent   = mgr->CreateContainer("ReducedEventDQ", AliReducedBaseEvent::Class(), AliAnalysisManager::kExchangeContainer, "reducedEvent");
    AliAnalysisDataContainer* cReducedTree    = 0x0;
    AliAnalysisDataContainer* cEventStatsInfo = 0x0;
    if(task->WriteTree()) {
      cReducedTree    = mgr->CreateContainer("dstTree",     TTree::Class(), AliAnalysisManager::kOutputContainer, "dstTree.root");
      cEventStatsInfo = mgr->CreateContainer("EventStats",  TList::Class(), AliAnalysisManager::kOutputContainer, "dstTree.root");
    }
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
  //
  // event cuts
  //
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  return eventCuts;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* CaloClusterFilter(Bool_t isAOD) {
  //
  // cluster cuts
  //
  AliDielectronClusterCuts* clusterCuts = new AliDielectronClusterCuts("clusters", "clusters");
  clusterCuts->SetCaloType(AliDielectronClusterCuts::kEMCal);
  clusterCuts->SetNCellsMinCut(2);
  clusterCuts->SetRejectExotics();
  clusterCuts->SetRequireTrackMatch(kTRUE);
  return clusterCuts;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* JPsiElectronFilter(Bool_t isAOD) {
  //
  // electron track cuts
  //
  AliDielectronCutGroup*  jpsiElectrons = new AliDielectronCutGroup("jpsiElectrons", "J/psi candidate electrons", AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts*   trackCuts     = new AliDielectronVarCuts("trackCuts", "track cuts");
  trackCuts->AddCut(AliDielectronVarManager::kImpactParXY,  -1.0, 1.0);
  trackCuts->AddCut(AliDielectronVarManager::kImpactParZ,   -3.0, 3.0);
  trackCuts->AddCut(AliDielectronVarManager::kEta,          -0.9, 0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,           1.0, 1.0e+30);
  trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,      70.0, 161.0);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaEle, -3.0, 4.0);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPro, -2000.0, -2.0, kTRUE);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPio, -2000.0, -2.0, kTRUE);
  jpsiElectrons->AddCut(trackCuts);
  AliDielectronTrackCuts* trackCuts1 = new AliDielectronTrackCuts("trackCuts1", "track cuts 1");
  trackCuts1->SetRequireITSRefit(kTRUE);
  trackCuts1->SetRequireTPCRefit(kTRUE);
  if (isAOD) trackCuts1->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  jpsiElectrons->AddCut(trackCuts1);
  return jpsiElectrons;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* JPsiElectronPreFilter(Bool_t isAOD) {
  //
  // electron prefilter track cuts
  //
  AliDielectronCutGroup*  jpsiElectrons = new AliDielectronCutGroup("jpsiElectronsPrefilterCut", "J/psi candidate electrons (prefilter cut)", AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts*   trackCuts     = new AliDielectronVarCuts("trackCuts", "track cuts");
  trackCuts->AddCut(AliDielectronVarManager::kImpactParXY,  -3.0, 3.0);
  trackCuts->AddCut(AliDielectronVarManager::kImpactParZ,   -10.0, 10.0);
  trackCuts->AddCut(AliDielectronVarManager::kEta,          -0.9, 0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,           0.5, 1.0e+30);
  trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,      50.0, 161.0);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaEle, -4.0, 4.0);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPro, -2000.0, 2.0, kTRUE);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPio, -2000.0, 2.0, kTRUE);
  jpsiElectrons->AddCut(trackCuts);
  AliDielectronTrackCuts* trackCuts1 = new AliDielectronTrackCuts("trackCuts1","track cuts 1");
  trackCuts1->SetRequireTPCRefit(kTRUE);
  if (isAOD) trackCuts1->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  jpsiElectrons->AddCut(trackCuts1);
  return jpsiElectrons;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* AssociatedHadronFilter(Bool_t isAOD, Int_t spec) {
  //
  // associated hadron track cuts
  // 0 - ITS + TPC refit
  // 1 - ITS + TPC refit, SPD any
  // 2 - ITS + TPC refit, no kinks
  // 3 - ITS + TPC refit, SPD any, no kinks
  //
  AliDielectronCutGroup*  assocHadr = new AliDielectronCutGroup(Form("assocHadron%d", spec), Form("Associated hadrons %d", spec), AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts*   trackCuts = new AliDielectronVarCuts("trackCuts", "track cuts");
  trackCuts->AddCut(AliDielectronVarManager::kImpactParXY,              -1.0, 1.0);
  trackCuts->AddCut(AliDielectronVarManager::kImpactParZ,               -3.0, 3.0);
  trackCuts->AddCut(AliDielectronVarManager::kEta,                      -0.9, 0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,                       0.15, 1.0e+30);
  trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,                  70.0, 161.0);
  if (spec>=2) trackCuts->AddCut(AliDielectronVarManager::kKinkIndex0,  -0.5, 0.5);
  assocHadr->AddCut(trackCuts);
  AliDielectronTrackCuts* trackCuts1 = new AliDielectronTrackCuts("trackCuts1", "track cuts 1");
  trackCuts1->SetRequireITSRefit(kTRUE);
  trackCuts1->SetRequireTPCRefit(kTRUE);
  if (spec==1 || spec==3) trackCuts1->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kAny);
  if (isAOD) trackCuts1->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  assocHadr->AddCut(trackCuts1);
  return assocHadr;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* AssociatedPionFilter(Bool_t isAOD, Int_t spec) {
  //
  // associated pion track cuts
  // 0 - ITS + TPC refit
  // 1 - ITS + TPC refit, SPD any
  // 2 - ITS + TPC refit, no kinks
  // 3 - ITS + TPC refit, SPD any, no kinks
  //
  AliDielectronCutGroup*  assocPi   = new AliDielectronCutGroup(Form("assocPion%d", spec), Form("Associated pions %d", spec), AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts*   trackCuts = new AliDielectronVarCuts("trackCuts", "track cuts");
  trackCuts->AddCut(AliDielectronVarManager::kImpactParXY,              -1.0, 1.0);
  trackCuts->AddCut(AliDielectronVarManager::kImpactParZ,               -3.0, 3.0);
  trackCuts->AddCut(AliDielectronVarManager::kEta,                      -0.9, 0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,                       0.15, 1.0e+30);
  trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,                  70.0, 161.0);
  if (spec>=2) trackCuts->AddCut(AliDielectronVarManager::kKinkIndex0,  -0.5, 0.5);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPio,             -3.0, 3.0);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaEle,             -2.0, 2000., kTRUE);
  assocPi->AddCut(trackCuts);
  AliDielectronTrackCuts* trackCuts1 = new AliDielectronTrackCuts("trackCuts1", "track cuts 1");
  trackCuts1->SetRequireITSRefit(kTRUE);
  trackCuts1->SetRequireTPCRefit(kTRUE);
  if (spec==1 || spec==3) trackCuts1->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kAny);
  if (isAOD) trackCuts1->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  assocPi->AddCut(trackCuts1);
  return assocPi;
}

//_______________________________________________________________________________________________________________
void AddMCSignals(AliAnalysisTaskReducedTreeMaker* task) {
  //
  // MC signals
  //

  // 0 = kJpsiInclusive
  AliSignalMC* jpsiInclusive = new AliSignalMC("JpsiInclusive", "",1,1);
  jpsiInclusive->SetPDGcode(0, 0, 443, kFALSE);
  task->AddMCsignal(jpsiInclusive, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 1 = kJpsiNonPrompt
  AliSignalMC* jpsiFromB = new AliSignalMC("JpsiNonPrompt","",1,2);
  jpsiFromB->SetPDGcode(0, 0, 443, kFALSE);
  jpsiFromB->SetPDGcode(0, 1, 503, kTRUE);
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
void SetInactiveBranches(AliAnalysisTaskReducedTreeMaker *task) {
  //
  // set inactive branches for treee
  //

  // event
  task->SetTreeInactiveBranch("fCentrality*");
  task->SetTreeInactiveBranch("fCentQuality");
  task->SetTreeInactiveBranch("fNV0candidates*");
  task->SetTreeInactiveBranch("fCandidates.*");
  task->SetTreeInactiveBranch("fEventNumberInFile");
  task->SetTreeInactiveBranch("fL0TriggerInputs");
  task->SetTreeInactiveBranch("fL1TriggerInputs");
  task->SetTreeInactiveBranch("fL2TriggerInputs");
  task->SetTreeInactiveBranch("fTimeStamp");
  task->SetTreeInactiveBranch("fEventType");
  task->SetTreeInactiveBranch("fMultiplicityEstimators*");
  task->SetTreeInactiveBranch("fMultiplicityEstimatorPercentiles*");
  task->SetTreeInactiveBranch("fIRIntClosestIntMap*");
  task->SetTreeInactiveBranch("fVtxTPC*");
  task->SetTreeInactiveBranch("fNVtxTPCContributors");
  task->SetTreeInactiveBranch("fVtxSPD*");
  task->SetTreeInactiveBranch("fNVtxSPDContributors");
  task->SetTreeInactiveBranch("fNPMDtracks");
  task->SetTreeInactiveBranch("fNTRDtracks");
  task->SetTreeInactiveBranch("fNTRDtracklets");
  task->SetTreeInactiveBranch("fSPDntrackletsEta*");
  task->SetTreeInactiveBranch("fVZEROMult*");
  task->SetTreeInactiveBranch("fZDCnEnergy*");
  task->SetTreeInactiveBranch("fZDCpEnergy*");
  task->SetTreeInactiveBranch("fZDCnTotalEnergy*");
  task->SetTreeInactiveBranch("fZDCpTotalEnergy*");
  task->SetTreeInactiveBranch("fT0amplitude*");
  task->SetTreeInactiveBranch("fT0zVertex");
  task->SetTreeInactiveBranch("fT0sattelite");
  task->SetTreeInactiveBranch("fFMD.*");
  task->SetTreeInactiveBranch("fEventPlane.*");
  task->SetTreeInactiveBranch("fTRDfired");
  task->SetTreeInactiveBranch("fBC");
  task->SetTreeInactiveBranch("fNtracksPerTrackingFlag*");
  task->SetTreeInactiveBranch("fT0TOF*");
  task->SetTreeInactiveBranch("fT0start");
  task->SetTreeInactiveBranch("fT0pileup");
  task->SetTreeInactiveBranch("fITSClusters*");
  task->SetTreeInactiveBranch("fSPDnSingle");

  // tracks
  task->SetTreeInactiveBranch("fTracks.fTPCPhi");
  task->SetTreeInactiveBranch("fTracks.fTPCPt");
  task->SetTreeInactiveBranch("fTracks.fTPCEta");
  task->SetTreeInactiveBranch("fTracks.fTPCDCA*");
  task->SetTreeInactiveBranch("fTracks.fTrackLength");
  task->SetTreeInactiveBranch("fTracks.fMassForTracking");
  task->SetTreeInactiveBranch("fTracks.fHelixCenter*");
  task->SetTreeInactiveBranch("fTracks.fHelixRadius");
  task->SetTreeInactiveBranch("fTracks.fITSsignal");
  task->SetTreeInactiveBranch("fTracks.fITSnSig*");
  task->SetTreeInactiveBranch("fTracks.fTPCActiveLength");
  task->SetTreeInactiveBranch("fTracks.fTPCGeomLength");
  task->SetTreeInactiveBranch("fTracks.fTPCdEdxInfoQmax*");
  task->SetTreeInactiveBranch("fTracks.fTPCdEdxInfoQtot*");
  task->SetTreeInactiveBranch("fTracks.fTRDntracklets*");
  task->SetTreeInactiveBranch("fTracks.fTRDpid*");
  task->SetTreeInactiveBranch("fTracks.fTRDpidLQ2D*");
  task->SetTreeInactiveBranch("fTracks.fTRDGTUtracklets");
  task->SetTreeInactiveBranch("fTracks.fTRDGTUlayermask");
  task->SetTreeInactiveBranch("fTracks.fTRDGTUpt");
  task->SetTreeInactiveBranch("fTracks.fTRDGTUsagitta");
  task->SetTreeInactiveBranch("fTracks.fTRDGTUPID");
  task->SetTreeInactiveBranch("fTracks.fTOFbeta");
  task->SetTreeInactiveBranch("fTracks.fTOFtime");
  task->SetTreeInactiveBranch("fTracks.fTOFdx");
  task->SetTreeInactiveBranch("fTracks.fTOFdz");
  task->SetTreeInactiveBranch("fTracks.fTOFmismatchProbab");
  task->SetTreeInactiveBranch("fTracks.fTOFchi2");
  task->SetTreeInactiveBranch("fTracks.fTOFnSig*");
  task->SetTreeInactiveBranch("fTracks.fTOFdeltaBC");
}
