void AddFMDTask();

//__________________________________________________________________________________________
AliAnalysisTask *AddTask_iarsene_dst(Int_t reducedEventType=-1, Bool_t writeTree=kTRUE, TString prod="LHC10h"){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_iarsene_dst", "No analysis manager found.");
    return 0;
  }

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  
  Int_t collSystem = 0;
  if(!prod.CompareTo("LHC10h")) collSystem = 1;    // minimum bias Pb-Pb
  if(!prod.CompareTo("LHC11h")) collSystem = 2;    // centrality triggered Pb-Pb
  if(!prod.CompareTo("LHC15o")) collSystem = 3;    // minimum bias Pb-Pb
  if(!prod.CompareTo("LHC10b") || !prod.CompareTo("LHC10c") || !prod.CompareTo("LHC10d") || 
     !prod.CompareTo("LHC10e") || !prod.CompareTo("LHC10f")) collSystem = 4;   // minimum bias pp
  if(!prod.CompareTo("LHC13b") || !prod.CompareTo("LHC13c")) collSystem = 5;   // minimum bias p-Pb   
  if(!prod.CompareTo("LHC16l")) collSystem = 6;    // minimum bias pp and Pb-Pb
  if(!prod.CompareTo("LHC16q") || !prod.CompareTo("LHC16r") || !prod.CompareTo("LHC16s") || !prod.CompareTo("LHC16t")) collSystem = 7;    // triggered p-Pb (2016 data)
  if(!prod.CompareTo("LHC17n")) collSystem = 8;
  
  if(!collSystem && !hasMC) {
     collSystem = 1;
     printf("WARNING: In AddTask_iarsene_dst(), no proper production name specified, or not supported! \n");
     printf("                 Using collSystem=1 (min bias 2010 Pb-Pb) as default \n");
  }
  
  //create task and add it to the manager
  AliAnalysisTaskReducedTreeMaker *task=new AliAnalysisTaskReducedTreeMaker("DSTTreeMaker", kTRUE);
  if(collSystem==1) task->SetTriggerMask(AliVEvent::kMB);
  if(collSystem==2)
    task->SetTriggerMask(AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral);
  if(collSystem==3) task->SetTriggerMask(AliVEvent::kINT7);
  if(collSystem==4) task->SetTriggerMask(AliVEvent::kMB);
  if(collSystem==5) task->SetTriggerMask(AliVEvent::kINT7);
  if(collSystem==6) task->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kHighMultV0);
  //if(collSystem==6) task->SetTriggerMask(AliVEvent::kINT7);
  //if(collSystem==6) task->SetTriggerMask(AliVEvent::kHighMultV0);
  //if(collSystem==6) task->SetTriggerMask(AliVEvent::kAny);
  if(collSystem==7) task->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kTRD);
  if(collSystem==8) task->SetTriggerMask(AliVEvent::kINT7);
  
  //if(collSystem==4)
  //  task->SetRejectPileup();
    
  task->UsePhysicsSelection(kTRUE);
  task->SetUseAnalysisUtils(kTRUE);
  
  task->SetFillV0Info(kFALSE);
  //task->SetFillGammaConversions(kFALSE);
  //task->SetFillK0s(kFALSE);
  //task->SetFillLambda(kFALSE);
  //task->SetFillALambda(kFALSE);
  //task->SetFillCaloClusterInfo(kTRUE);
  //task->SetFillDielectronInfo(kFALSE);
  //task->SetFillFriendInfo(kFALSE);
  
  task->SetEventFilter(CreateEventFilter(isAOD));
  //task->SetTrackFilter(CreateGlobalTrackFilter(isAOD));
  CreateTrackFilters(task, isAOD);
  
  //  task->SetFlowTrackFilter(CreateFlowTrackFilter(isAOD));
  //task->SetK0sPionCuts(CreateK0sPionCuts(isAOD));
  //task->SetLambdaProtonCuts(CreateLambdaProtonCuts(isAOD));
  //task->SetLambdaPionCuts(CreateLambdaPionCuts(isAOD));
  //task->SetGammaElectronCuts(CreateGammaConvElectronCuts(isAOD));
  //task->SetK0sMassRange(0.44,0.55);
  //task->SetLambdaMassRange(1.090,1.14);
  //task->SetGammaConvMassRange(0.0,0.1);
  task->SetV0OpenCuts(CreateV0OpenCuts(AliESDv0KineCuts::kPurity, AliESDv0KineCuts::kPbPb));
  task->SetV0StrongCuts(CreateV0StrongCuts(AliESDv0KineCuts::kPurity, AliESDv0KineCuts::kPbPb));
  //task->SetFillFMDInfo(); AddFMDTask();
  if(hasMC) {
     task->SetFillMCInfo(kTRUE);
     AddMCSignals(task);
  }
  //task->SetFillEventPlaneInfo(kTRUE);
  
  task->SetWriteEventsWithNoSelectedTracks(kFALSE,0.01,2);
  task->SetWriteEventsWithNoSelectedTracksAndNoSelectedAssociatedTracks(kFALSE);

  //task->SetTreeWritingOption(AliAnalysisTaskReducedTreeMaker::kFullEventsWithFullTracks);
  //task->SetTreeWritingOption(AliAnalysisTaskReducedTreeMaker::kBaseEventsWithBaseTracks);
  if(reducedEventType!=-1)
    task->SetTreeWritingOption(reducedEventType);
  task->SetWriteTree(writeTree);
  //task->SetWriteEventsWithNoSelectedTracks(kFALSE);
  
  SetInactiveBranches(task);
  
  mgr->AddTask(task);
    
  AliAnalysisDataContainer *cReducedTree = 0x0;
  AliAnalysisDataContainer *cEventStatsInfo = 0x0;
  AliAnalysisDataContainer *cTrackStatsInfo = 0x0;
  AliAnalysisDataContainer *cMCStatsInfo = 0x0;
  if(task->WriteTree()) {
    cReducedTree = mgr->CreateContainer("dstTree", TTree::Class(),
                                          AliAnalysisManager::kOutputContainer, "dstTree.root");
    cEventStatsInfo = mgr->CreateContainer("EventStats", TH2I::Class(),
       AliAnalysisManager::kOutputContainer, "dstTree.root");
    cTrackStatsInfo = mgr->CreateContainer("TrackStats", TH2I::Class(),
                                           AliAnalysisManager::kOutputContainer, "dstTree.root");
    cMCStatsInfo = mgr->CreateContainer("MCStats", TH2I::Class(),
                                           AliAnalysisManager::kOutputContainer, "dstTree.root");
  }
  
  AliAnalysisDataContainer *cReducedEvent =
  mgr->CreateContainer("ReducedEventDQ",
                         AliReducedBaseEvent::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         "reducedEvent");
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cReducedEvent);
  if(task->WriteTree()) {
     mgr->ConnectOutput(task, 2, cReducedTree);
     mgr->ConnectOutput(task, 3, cEventStatsInfo);
     mgr->ConnectOutput(task, 4, cTrackStatsInfo);
     if(hasMC) mgr->ConnectOutput(task, 5, cMCStatsInfo);
  }
  
  return task;
}

//_______________________________________________________________________________________________
void AddMCSignals(AliAnalysisTaskReducedTreeMaker* task) {
   //
   // Add the MC signals to be filtered
   //
   AliSignalMC* jpsiInclusive=new AliSignalMC("JpsiInclusive", "",1,1);
   jpsiInclusive->SetPDGcode(0, 0, 443, kFALSE);
   task->AddMCsignal(jpsiInclusive, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   AliSignalMC* jpsiFromB=new AliSignalMC("JpsiFromB","",1,2);
   jpsiFromB->SetPDGcode(0, 0, 443, kFALSE);
   jpsiFromB->SetPDGcode(0, 1, 500, kTRUE);
   task->AddMCsignal(jpsiFromB, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   AliSignalMC* jpsiPrompt=new AliSignalMC("JpsiNotFromB","",1,2);
   jpsiPrompt->SetPDGcode(0, 0, 443, kFALSE);
   jpsiPrompt->SetPDGcode(0, 1, 500, kTRUE, kTRUE);
   task->AddMCsignal(jpsiPrompt, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   AliSignalMC* jpsiInclusiveRadiative=new AliSignalMC("JpsiInclusiveRadiative","",1,1);
   jpsiInclusiveRadiative->SetPDGcode(0, 0, 443, kFALSE);
   jpsiInclusiveRadiative->SetSourceBit(0, 0, AliSignalMC::kRadiativeDecay);
   task->AddMCsignal(jpsiInclusiveRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   AliSignalMC* jpsiInclusiveNonRadiative=new AliSignalMC("JpsiInclusiveNonRadiative","",1,1);
   jpsiInclusiveNonRadiative->SetPDGcode(0, 0, 443, kFALSE);
   jpsiInclusiveNonRadiative->SetSourceBit(0, 0, AliSignalMC::kRadiativeDecay, kTRUE);
   task->AddMCsignal(jpsiInclusiveNonRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   AliSignalMC* jpsiFromBRadiative=new AliSignalMC("JpsiFromBRadiative","",1,2);
   jpsiFromBRadiative->SetPDGcode(0, 0, 443, kFALSE);
   jpsiFromBRadiative->SetPDGcode(0, 1, 500, kTRUE);
   jpsiFromBRadiative->SetSourceBit(0, 0, AliSignalMC::kRadiativeDecay, kFALSE);
   task->AddMCsignal(jpsiFromBRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   AliSignalMC* jpsiFromBNonRadiative=new AliSignalMC("JpsiFromBNonRadiative","",1,2);
   jpsiFromBNonRadiative->SetPDGcode(0, 0, 443, kFALSE);
   jpsiFromBNonRadiative->SetPDGcode(0, 1, 500, kTRUE);
   jpsiFromBNonRadiative->SetSourceBit(0, 0, AliSignalMC::kRadiativeDecay, kTRUE);
   task->AddMCsignal(jpsiFromBNonRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   
   AliSignalMC* electronFromJpsi=new AliSignalMC(1,2);
   electronFromJpsi->SetPDGcode(0, 0, 11, kTRUE);
   electronFromJpsi->SetPDGcode(0, 1, 443);
   //task->AddMCsignal(electronFromJpsi, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   AliSignalMC* electronFromJpsiPrompt=new AliSignalMC(1,3);
   electronFromJpsiPrompt->SetPDGcode(0, 0, 11, kTRUE);
   electronFromJpsiPrompt->SetPDGcode(0, 1, 443);
   electronFromJpsiPrompt->SetPDGcode(0, 2, 500, kTRUE, kTRUE);    //
   //task->AddMCsignal(electronFromJpsiPrompt, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   AliSignalMC* electronFromJpsiNonPrompt=new AliSignalMC(1,3);
   electronFromJpsiNonPrompt->SetPDGcode(0, 0, 11, kTRUE);
   electronFromJpsiNonPrompt->SetPDGcode(0, 1, 443);
   electronFromJpsiNonPrompt->SetPDGcode(0, 2, 500, kTRUE);      //
   //task->AddMCsignal(electronFromJpsiNonPrompt, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   AliSignalMC* photonFromJpsiDecay=new AliSignalMC(1,2);
   photonFromJpsiDecay->SetPDGcode(0, 0, 22);
   photonFromJpsiDecay->SetPDGcode(0, 1, 443);
   //task->AddMCsignal(photonFromJpsiDecay, AliAnalysisTaskReducedTreeMaker::kFullTrack);
}

//_______________________________________________________________________________________________
void SetInactiveBranches(AliAnalysisTaskReducedTreeMaker* task) {
   //
   // Set inactive branches
   //
   // You can define the active branches of the tree, all other branches will be set to inactive
   //task->SetTreeActiveBranch("Event");
   //task->SetTreeActiveBranch("fRunNo");
   // Alternatively, you can define the inactive branches of the tree, all other branches will be set to active

   
   //task->SetTreeInactiveBranch("fEventTag");
   //task->SetTreeInactiveBranch("fRunNo");
   //task->SetTreeInactiveBranch("fVtx*");
   //task->SetTreeInactiveBranch("fNVtxContributors");
   task->SetTreeInactiveBranch("fCentrality*");
   task->SetTreeInactiveBranch("fCentQuality");
   //task->SetTreeInactiveBranch("fNtracks*");
   task->SetTreeInactiveBranch("fNV0candidates*");
   //task->SetTreeInactiveBranch("fTracks.*");
   task->SetTreeInactiveBranch("fCandidates.*");
   task->SetTreeInactiveBranch("fEventNumberInFile");
   task->SetTreeInactiveBranch("fL0TriggerInputs");
   task->SetTreeInactiveBranch("fL1TriggerInputs");
   task->SetTreeInactiveBranch("fL2TriggerInputs");
   task->SetTreeInactiveBranch("fBC");
   task->SetTreeInactiveBranch("fTimeStamp");
   task->SetTreeInactiveBranch("fEventType");
   //task->SetTreeInactiveBranch("fTriggerMask");
   task->SetTreeInactiveBranch("fMultiplicityEstimators*");
   task->SetTreeInactiveBranch("fMultiplicityEstimatorPercentiles*");
   //task->SetTreeInactiveBranch("fIsPhysicsSelection");
   //task->SetTreeInactiveBranch("fIsSPDPileup");
   //task->SetTreeInactiveBranch("fIsSPDPileupMultBins");
   task->SetTreeInactiveBranch("fIRIntClosestIntMap*");
   task->SetTreeInactiveBranch("fVtxTPC*");
   task->SetTreeInactiveBranch("fNVtxTPCContributors");
   task->SetTreeInactiveBranch("fVtxSPD*");
   task->SetTreeInactiveBranch("fNVtxSPDContributors");
   task->SetTreeInactiveBranch("fNpileupSPD");
   task->SetTreeInactiveBranch("fNpileupTracks");
   //task->SetTreeInactiveBranch("fNTPCclusters");
   task->SetTreeInactiveBranch("fNPMDtracks");
   task->SetTreeInactiveBranch("fNTRDtracks");
   task->SetTreeInactiveBranch("fNTRDtracklets");
   //task->SetTreeInactiveBranch("fSPDntracklets");
   task->SetTreeInactiveBranch("fSPDntrackletsEta*");
   //task->SetTreeInactiveBranch("fSPDFiredChips*");
   task->SetTreeInactiveBranch("fITSClusters*");
   //task->SetTreeInactiveBranch("fSPDnSingle");
   task->SetTreeInactiveBranch("fNtracksPerTrackingFlag*");
   task->SetTreeInactiveBranch("fVZEROMult*");
   //task->SetTreeInactiveBranch("fVZEROTotalMult*");
   task->SetTreeInactiveBranch("fZDCnEnergy*");
   task->SetTreeInactiveBranch("fZDCpEnergy*");
   task->SetTreeInactiveBranch("fZDCnTotalEnergy*");
   task->SetTreeInactiveBranch("fZDCpTotalEnergy*");
   task->SetTreeInactiveBranch("fT0amplitude*");
   task->SetTreeInactiveBranch("fT0TOF*");
   task->SetTreeInactiveBranch("fT0TOFbest*");
   task->SetTreeInactiveBranch("fT0zVertex");
   task->SetTreeInactiveBranch("fT0start");
   task->SetTreeInactiveBranch("fT0pileup");
   task->SetTreeInactiveBranch("fT0sattelite");
   task->SetTreeInactiveBranch("fNCaloClusters");
   task->SetTreeInactiveBranch("fCaloClusters.*");
   task->SetTreeInactiveBranch("fFMD.*");
   task->SetTreeInactiveBranch("fEventPlane.*");
   
   //task->SetTreeInactiveBranch("fTracks2.*");
   task->SetTreeInactiveBranch("fTracks2.fFlags");
   task->SetTreeInactiveBranch("fTracks2.fMCFlags");
   //task->SetTreeInactiveBranch("fTracks.fP*");
   //task->SetTreeInactiveBranch("fTracks.fIsCartesian");
   //task->SetTreeInactiveBranch("fTracks.fCharge");
   task->SetTreeInactiveBranch("fTracks.fFlags");
   task->SetTreeInactiveBranch("fTracks.fMCFlags");
   //task->SetTreeInactiveBranch("fTracks.fQualityFlags");
   //task->SetTreeInactiveBranch("fTracks.fTrackId");
   //task->SetTreeInactiveBranch("fTracks.fStatus");
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
   //task->SetTreeInactiveBranch("fTracks.fITSchi2");
   //task->SetTreeInactiveBranch("fTracks.fTPCNcls");
   //task->SetTreeInactiveBranch("fTracks.fTPCCrossedRows");
   //task->SetTreeInactiveBranch("fTracks.fTPCNclsF");
   //task->SetTreeInactiveBranch("fTracks.fTPCNclsShared");
   //task->SetTreeInactiveBranch("fTracks.fTPCClusterMap");
   //task->SetTreeInactiveBranch("fTracks.fTPCsignal");
   //task->SetTreeInactiveBranch("fTracks.fTPCsignalN");
   //task->SetTreeInactiveBranch("fTracks.fTPCnSig*");
   //task->SetTreeInactiveBranch("fTracks.fTPCchi2");
   task->SetTreeInactiveBranch("fTracks.fTPCActiveLength");
   task->SetTreeInactiveBranch("fTracks.fTPCGeomLength");
   //task->SetTreeInactiveBranch("fTracks.fTOFbeta");
   task->SetTreeInactiveBranch("fTracks.fTOFtime");
   task->SetTreeInactiveBranch("fTracks.fTOFdx");
   task->SetTreeInactiveBranch("fTracks.fTOFdz");
   task->SetTreeInactiveBranch("fTracks.fTOFmismatchProbab");
   task->SetTreeInactiveBranch("fTracks.fTOFchi2");
   task->SetTreeInactiveBranch("fTracks.fTOFnSig*");
   task->SetTreeInactiveBranch("fTracks.fTOFdeltaBC");
   task->SetTreeInactiveBranch("fTracks.fTRDntracklets*");
   task->SetTreeInactiveBranch("fTracks.fTRDpid*");
   task->SetTreeInactiveBranch("fTracks.fTRDpidLQ2D*");
   task->SetTreeInactiveBranch("fTracks.fCaloClusterId");
   task->SetTreeInactiveBranch("fTracks.fMCMom*");
   task->SetTreeInactiveBranch("fTracks.fMCFreezeout*");
   task->SetTreeInactiveBranch("fTracks.fMCLabels*");
   task->SetTreeInactiveBranch("fTracks.fMCPdg*");
   task->SetTreeInactiveBranch("fTracks.fMCGeneratorIndex");
  
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
  trackCuts->AddCut(AliDielectronVarManager::kP,1.0,1.0e+30);
  trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,70.0,161.0);
  //  trackCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,0.1,4.0);
  
  //trackCuts->AddCut(AliDielectronVarManager::kP,1.0,1.0e+30);
  //trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,70.0,161.0);
  cuts->AddCut(trackCuts);
  AliDielectronTrackCuts* trackCuts2 = new AliDielectronTrackCuts("trackCuts2","track cuts");
  //trackCuts2->SetRequireITSRefit(kTRUE);
  trackCuts2->SetRequireTPCRefit(kTRUE);
  if(isAOD) trackCuts2->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  cuts->AddCut(trackCuts2);
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
  //electronPid->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-2.0, 2.0, 0.0, 0.0, kTRUE, AliDielectronPID::kRequire); // TPC 3-sigma inclusion for proton
  //electronPid->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,-3.0, 3.0, 0.0, 0.0, kTRUE, AliDielectronPID::kRequire); // TPC 3-sigma inclusion for kaon
  //electronPid->AddCut(AliDielectronPID::kTOF,AliPID::kProton,  -3.0, 3.0, -2.0, 2.0, kTRUE, AliDielectronPID::kIfAvailable, AliDielectronVarManager::kTPCnSigmaPro); // TPC exclusion for proton   
  cuts->AddCut(electronPid);
  
  return cuts;
}


//______________________________________________________________________________________
void CreateTrackFilters(AliAnalysisTaskReducedTreeMaker* task, Bool_t isAOD) {
   //
   // add track filters
   //
   AliDielectronCutGroup* jpsiElectrons = new AliDielectronCutGroup("jpsiElectrons","J/psi candidate electrons",AliDielectronCutGroup::kCompAND);
   AliDielectronVarCuts *trackCuts = new AliDielectronVarCuts("trackCuts","track cuts");
   trackCuts->AddCut(AliDielectronVarManager::kImpactParXY,-1.0,1.0);
   trackCuts->AddCut(AliDielectronVarManager::kImpactParZ,-3.0,3.0);
   trackCuts->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
   trackCuts->AddCut(AliDielectronVarManager::kPt,1.0,1.0e+30);
   trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,70.0,161.0);   
   trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaEle, -3.0, 3.0);   
   trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPro, -2000.0, 3.0, kTRUE);
   trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPio, -2000.0, 2.0, kTRUE);
   jpsiElectrons->AddCut(trackCuts);
   AliDielectronTrackCuts* trackCuts2 = new AliDielectronTrackCuts("trackCuts2","track cuts");
   trackCuts2->SetRequireITSRefit(kTRUE);
   trackCuts2->SetRequireTPCRefit(kTRUE);
   if(isAOD) trackCuts2->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
   jpsiElectrons->AddCut(trackCuts2);
   task->AddTrackFilter(jpsiElectrons, kFALSE);
   
   AliDielectronCutGroup* assocHadr = new AliDielectronCutGroup("assocHadr","Associated hadrons",AliDielectronCutGroup::kCompAND);
   AliDielectronVarCuts *trackCuts3 = new AliDielectronVarCuts("trackCuts3","track cuts3");
   trackCuts3->AddCut(AliDielectronVarManager::kImpactParXY,-1.0,1.0);
   trackCuts3->AddCut(AliDielectronVarManager::kImpactParZ,-3.0,3.0);
   trackCuts3->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
   trackCuts3->AddCut(AliDielectronVarManager::kPt,0.5,1.0e+30);
   trackCuts3->AddCut(AliDielectronVarManager::kNclsTPC,70.0,161.0);   
   assocHadr->AddCut(trackCuts3);
   AliDielectronTrackCuts* trackCuts4 = new AliDielectronTrackCuts("trackCuts4","track cuts4");
   trackCuts4->SetRequireITSRefit(kTRUE);
   trackCuts4->SetRequireTPCRefit(kTRUE);
   if(isAOD) trackCuts4->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
   assocHadr->AddCut(trackCuts4);
   task->AddTrackFilter(assocHadr, kTRUE);
   
   AliDielectronCutGroup* assocHadr2 = new AliDielectronCutGroup("assocHadr2","Associated hadrons, no ITSrefit",AliDielectronCutGroup::kCompAND);
   assocHadr2->AddCut(trackCuts3);
   AliDielectronTrackCuts* trackCuts5 = new AliDielectronTrackCuts("trackCuts5","track cuts5");
   trackCuts5->SetRequireTPCRefit(kTRUE);
   if(isAOD) trackCuts5->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
   assocHadr2->AddCut(trackCuts5);
   task->AddTrackFilter(assocHadr2, kTRUE);
   
   AliDielectronCutGroup* assocPion = new AliDielectronCutGroup("assocPion","Associated pions",AliDielectronCutGroup::kCompAND);
   assocPion->AddCut(trackCuts3);
   assocPion->AddCut(trackCuts4);
   AliDielectronVarCuts *pionPID = new AliDielectronVarCuts("pionPID","pion PID");
   pionPID->AddCut(AliDielectronVarManager::kTPCnSigmaPio, -3.0, 3.0);
   assocPion->AddCut(pionPID);
   task->AddTrackFilter(assocPion, kTRUE);
   
   AliDielectronCutGroup* assocKaon = new AliDielectronCutGroup("assocKaon","Associated kaons",AliDielectronCutGroup::kCompAND);
   assocKaon->AddCut(trackCuts3);
   assocKaon->AddCut(trackCuts4);
   AliDielectronVarCuts *kaonPID = new AliDielectronVarCuts("kaonPID","kaon PID");
   kaonPID->AddCut(AliDielectronVarManager::kTPCnSigmaKao, -3.0, 3.0);
   assocKaon->AddCut(kaonPID);
   task->AddTrackFilter(assocKaon, kTRUE);
   
   AliDielectronCutGroup* assocProton = new AliDielectronCutGroup("assocProton","Associated protons",AliDielectronCutGroup::kCompAND);
   assocProton->AddCut(trackCuts3);
   assocProton->AddCut(trackCuts4);
   AliDielectronVarCuts *protonPID = new AliDielectronVarCuts("protonPID","proton PID");
   protonPID->AddCut(AliDielectronVarManager::kTPCnSigmaPro, -3.0, 3.0);
   assocProton->AddCut(protonPID);
   task->AddTrackFilter(assocProton, kTRUE);
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
