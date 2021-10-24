//===============================================================================================================
// addtask to create trees for single track hadron efficiency in pp 13TeV (last updated: 10/11/2019)
//===============================================================================================================

AliAnalysisCuts* EventFilter(Bool_t isAOD);
AliAnalysisCuts* AssociatedHadronFilter(Bool_t isAOD, Int_t spec);
AliAnalysisCuts* AssociatedPionFilter(Bool_t isAOD, Int_t spec);
void AddMCSignals(AliAnalysisTaskReducedTreeMaker* task);
void SetActiveBranches(AliAnalysisTaskReducedTreeMaker* task);

//_______________________________________________________________________________________________________________
AliAnalysisTask* AddTask_laltenka_dst_assocHadronEff(Int_t reducedEventType=-1, Bool_t writeTree=kTRUE, TString  prod="") {

  //get current analysis manager
  //-----------------------------------------------------------------------------------------------------------
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) { Error("AddTask_laltenka_dst_assocHadronEff", "No analysis manager found."); return 0; }
  
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
  printf("AddTask_laltenka_dst_assocHadronEff() trigger choice set to %d (%s)\n", triggerChoice, prod.Data());
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
    printf("WARNING: In AddTask_laltenka_dst_assocHadronEff(), no trigger specified, or not supported! Using kINT7.\n");
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
  task->SetFillV0Info(kFALSE);
  task->SetFillGammaConversions(kFALSE);
  task->SetFillK0s(kFALSE);
  task->SetFillLambda(kFALSE);
  task->SetFillALambda(kFALSE);
  task->SetFillCaloClusterInfo(kFALSE);
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

  // track selection
  //-----------------------------------------------------------------------------------------------------------
  // associated unidentified hadrons
  task->AddTrackFilter(AssociatedHadronFilter(isAOD, 0),  kFALSE, 1e6); // base track info, fQualityFlag bit 32 - unid. hadron, ITS + TPC refit
  task->AddTrackFilter(AssociatedHadronFilter(isAOD, 1),  kFALSE, 1e6); // base track info, fQualityFlag bit 33 - unid. hadron, ITS + TPC refit, SPD any
  task->AddTrackFilter(AssociatedHadronFilter(isAOD, 2),  kFALSE, 1e6); // base track info, fQualityFlag bit 34 - unid. hadron, ITS + TPC refit, no kinks
  task->AddTrackFilter(AssociatedHadronFilter(isAOD, 3),  kFALSE, 1e6); // base track info, fQualityFlag bit 35 - unid. hadron, ITS + TPC refit, SPD any, no kinks
  // associated pions
  task->AddTrackFilter(AssociatedPionFilter(isAOD, 0),    kFALSE, 1e6); // base track info, fQualityFlag bit 36 - pions, ITS + TPC refit
  task->AddTrackFilter(AssociatedPionFilter(isAOD, 1),    kFALSE, 1e6); // base track info, fQualityFlag bit 37 - pions, ITS + TPC refit, SPD any
  task->AddTrackFilter(AssociatedPionFilter(isAOD, 2),    kFALSE, 1e6); // base track info, fQualityFlag bit 38 - pions, ITS + TPC refit, no kinks
  task->AddTrackFilter(AssociatedPionFilter(isAOD, 3),    kFALSE, 1e6); // base track info, fQualityFlag bit 39 - pions, ITS + TPC refit, SPD any, no kinks

  // active/inactive branches of tree
  //-----------------------------------------------------------------------------------------------------------
  SetActiveBranches(task);
  
  // set writing options
  //  AliAnalysisTaskReducedTreeMaker::kBaseEventsWithBaseTracks = 0
  //  AliAnalysisTaskReducedTreeMaker::kBaseEventsWithFullTracks = 1
  //  AliAnalysisTaskReducedTreeMaker::kFullEventsWithBaseTracks = 2
  //  AliAnalysisTaskReducedTreeMaker::kFullEventsWithFullTracks = 3
  //-----------------------------------------------------------------------------------------------------------
  if(reducedEventType!=-1) task->SetTreeWritingOption(reducedEventType);
  task->SetWriteTree(writeTree);
  task->SetWriteUnbiasedEvents(0.05);

  // add task to manager
  //-----------------------------------------------------------------------------------------------------------
  mgr->AddTask(task);

  // connect output containers
  //-----------------------------------------------------------------------------------------------------------
  AliAnalysisDataContainer* cReducedEvent   = mgr->CreateContainer("ReducedEventDQ", AliReducedBaseEvent::Class(), AliAnalysisManager::kExchangeContainer, "reducedEvent");
  AliAnalysisDataContainer* cReducedTree    = 0x0;
  AliAnalysisDataContainer* cEventStatsInfo = 0x0;
  if(task->WriteTree()) {
    cReducedTree    = mgr->CreateContainer("dstTree",       TTree::Class(), AliAnalysisManager::kOutputContainer, "dstTree.root");
    cEventStatsInfo = mgr->CreateContainer("EventStats",    TList::Class(), AliAnalysisManager::kOutputContainer, "dstTree.root");
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
  AliDielectronEventCuts* eventCuts = new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  return eventCuts;
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
  // 1 = kPion
  AliSignalMC* pionPrim = new AliSignalMC("pionPrim", "", 1, 1);
  pionPrim->SetPDGcode(0, 0, 211, kTRUE);
  pionPrim->SetSourceBit(0, 0, AliSignalMC::kPhysicalPrimary);
  task->AddMCsignal(pionPrim, AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker::kBaseTrack);  //AliAnalysisTaskReducedTreeMaker::kFullTrack););

  AliSignalMC* pionSec = new AliSignalMC("pionSec", "", 1, 1);
  pionSec->SetPDGcode(0, 0, 211, kTRUE);
  pionSec->SetSourceBit(0, 0, AliSignalMC::kSecondaryFromWeakDecay);
  task->AddMCsignal(pionSec, AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker::kBaseTrack);  //AliAnalysisTaskReducedTreeMaker::kFullTrack););

  AliSignalMC* pionSecMat = new AliSignalMC("pionSecMat", "", 1, 1);
  pionSecMat->SetPDGcode(0, 0, 211, kTRUE);
  pionSecMat->SetSourceBit(0, 0, AliSignalMC::kSecondaryFromMaterial);
  task->AddMCsignal(pionSecMat, AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker::kBaseTrack);  //AliAnalysisTaskReducedTreeMaker::kFullTrack););

  // 2 = kKaon
  AliSignalMC* kaonPrim = new AliSignalMC("kaonPrim", "", 1, 1);
  kaonPrim->SetPDGcode(0, 0, 321, kTRUE);
  kaonPrim->SetSourceBit(0, 0, AliSignalMC::kPhysicalPrimary);
  task->AddMCsignal(kaonPrim, AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker::kBaseTrack);  //AliAnalysisTaskReducedTreeMaker::kFullTrack););

  AliSignalMC* kaonSec = new AliSignalMC("kaonSec", "", 1, 1);
  kaonSec->SetPDGcode(0, 0, 321, kTRUE);
  kaonSec->SetSourceBit(0, 0, AliSignalMC::kSecondaryFromWeakDecay);
  task->AddMCsignal(kaonSec, AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker::kBaseTrack);  //AliAnalysisTaskReducedTreeMaker::kFullTrack););

  AliSignalMC* kaonSecMat = new AliSignalMC("kaonSecMat", "", 1, 1);
  kaonSecMat->SetPDGcode(0, 0, 321, kTRUE);
  kaonSecMat->SetSourceBit(0, 0, AliSignalMC::kSecondaryFromMaterial);
  task->AddMCsignal(kaonSecMat, AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker::kBaseTrack);  //AliAnalysisTaskReducedTreeMaker::kFullTrack););

  // 3 = kProton
  AliSignalMC* protonPrim = new AliSignalMC("protonPrim", "", 1, 1);
  protonPrim->SetPDGcode(0, 0, 2212, kTRUE);
  protonPrim->SetSourceBit(0, 0, AliSignalMC::kPhysicalPrimary);
  task->AddMCsignal(protonPrim, AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker::kBaseTrack);  //AliAnalysisTaskReducedTreeMaker::kFullTrack););

  AliSignalMC* protonSec = new AliSignalMC("protonSec", "", 1, 1);
  protonSec->SetPDGcode(0, 0, 2212, kTRUE);
  protonSec->SetSourceBit(0, 0, AliSignalMC::kSecondaryFromWeakDecay);
  task->AddMCsignal(protonSec, AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker::kBaseTrack);  //AliAnalysisTaskReducedTreeMaker::kFullTrack););

  AliSignalMC* protonSecMat = new AliSignalMC("protonSecMat", "", 1, 1);
  protonSecMat->SetPDGcode(0, 0, 2212, kTRUE);
  protonSecMat->SetSourceBit(0, 0, AliSignalMC::kSecondaryFromMaterial);
  task->AddMCsignal(protonSecMat, AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker::kBaseTrack);  //AliAnalysisTaskReducedTreeMaker::kFullTrack););

  // 4 = kElectron
  AliSignalMC* electronPrim = new AliSignalMC("electronPrim", "", 1, 1);
  electronPrim->SetPDGcode(0, 0, 11, kTRUE);
  electronPrim->SetSourceBit(0, 0, AliSignalMC::kPhysicalPrimary);
  task->AddMCsignal(electronPrim, AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker::kBaseTrack);  //AliAnalysisTaskReducedTreeMaker::kFullTrack););

  AliSignalMC* electronSec = new AliSignalMC("electronSec", "", 1, 1);
  electronSec->SetPDGcode(0, 0, 11, kTRUE);
  electronSec->SetSourceBit(0, 0, AliSignalMC::kSecondaryFromWeakDecay);
  task->AddMCsignal(electronSec, AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker::kBaseTrack);  //AliAnalysisTaskReducedTreeMaker::kFullTrack););

  AliSignalMC* electronSecMat = new AliSignalMC("electronSecMat", "", 1, 1);
  electronSecMat->SetPDGcode(0, 0, 11, kTRUE);
  electronSecMat->SetSourceBit(0, 0, AliSignalMC::kSecondaryFromMaterial);
  task->AddMCsignal(electronSecMat, AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker::kBaseTrack);  //AliAnalysisTaskReducedTreeMaker::kFullTrack););

  // 5 = kMuon
  AliSignalMC* muonPrim = new AliSignalMC("muonPrim", "", 1, 1);
  muonPrim->SetPDGcode(0, 0, 13, kTRUE);
  muonPrim->SetSourceBit(0, 0, AliSignalMC::kPhysicalPrimary);
  task->AddMCsignal(muonPrim, AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker::kBaseTrack);  //AliAnalysisTaskReducedTreeMaker::kFullTrack););

  AliSignalMC* muonSec = new AliSignalMC("muonSec", "", 1, 1);
  muonSec->SetPDGcode(0, 0, 13, kTRUE);
  muonSec->SetSourceBit(0, 0, AliSignalMC::kSecondaryFromWeakDecay);
  task->AddMCsignal(muonSec, AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker::kBaseTrack);  //AliAnalysisTaskReducedTreeMaker::kFullTrack););

  AliSignalMC* muonSecMat = new AliSignalMC("muonSecMat", "", 1, 1);
  muonSecMat->SetPDGcode(0, 0, 13, kTRUE);
  muonSecMat->SetSourceBit(0, 0, AliSignalMC::kSecondaryFromMaterial);
  task->AddMCsignal(muonSecMat, AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker::kBaseTrack);  //AliAnalysisTaskReducedTreeMaker::kFullTrack););
}

//_______________________________________________________________________________________________________________
void SetActiveBranches(AliAnalysisTaskReducedTreeMaker* task) {
  //
  // set active branches for treee
  //
  
  // event
  task->SetTreeActiveBranch("fIsPhysicsSelection");
  task->SetTreeActiveBranch("fTriggerMask");
  task->SetTreeActiveBranch("fTriggerClass");
  task->SetTreeActiveBranch("fEventTag");
  task->SetTreeActiveBranch("fRunNo");
  task->SetTreeActiveBranch("fVtx*");
  task->SetTreeInactiveBranch("fVtxMC*");

  // track
  task->SetTreeActiveBranch("fTracks.fUniqueID");
  task->SetTreeActiveBranch("fTracks.fBits");
  task->SetTreeActiveBranch("fTracks.fTrackId");
  task->SetTreeActiveBranch("fTracks.fFlags");
  task->SetTreeActiveBranch("fTracks.fMCFlags");
  task->SetTreeActiveBranch("fTracks.fIsMCTruth");
  task->SetTreeActiveBranch("fTracks.fQualityFlags");
  task->SetTreeActiveBranch("fTracks.fStatus");
  task->SetTreeActiveBranch("fTracks.fP*");
  task->SetTreeActiveBranch("fTracks.fMCMom*");
  task->SetTreeActiveBranch("fTracks.fMCLabels*");
  task->SetTreeActiveBranch("fTracks.fMCPdg*");
  task->SetTreeActiveBranch("fTracks.fIsCartesian");
  task->SetTreeActiveBranch("fTracks.fCharge");
  task->SetTreeActiveBranch("fTracks.fDCA*");
}
