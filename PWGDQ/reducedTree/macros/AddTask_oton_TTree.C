//void AddFMDTask(); //is this actually necessary???

//__________________________________________________________________________________________
AliAnalysisTask *AddTask_oton_TTree(Int_t reducedEventType=-1, Bool_t writeTree=kTRUE){

  // Get the current analysis manager
  //---------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskTree", "No analysis manager found.");
    return 0;
  }

  // Do we have an MC handler?
  //--------------------------
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  
  //create task and add it to the manager
  //-------------------------------------
  AliAnalysisTaskReducedTreeMaker *task=new AliAnalysisTaskReducedTreeMaker("DSTTreeMaker", kTRUE);
  task->SetTriggerMask(AliVEvent::kINT7);
  
  // Pile up, Physics Selection, Analysis Utils
  //-------------------------------------------
  //if(trainConfig=="pp") task->SetRejectPileup();
  //task->UsePhysicsSelection(kFALSE); // No phys sel // is this necessary???
  task->SetUseAnalysisUtils(kFALSE); // moved to false 16/10/2016 // ???
  
  // Toggle on/off information branches
  //----------------------------------------
  task->SetFillV0Info(kFALSE);
  task->SetFillGammaConversions(kFALSE);
  task->SetFillK0s(kFALSE);
  task->SetFillLambda(kFALSE);
  task->SetFillALambda(kFALSE);
  task->SetFillCaloClusterInfo(kFALSE);
  //task->SetFillDielectronInfo(kFALSE);
  //task->SetFillFriendInfo(kFALSE);
  //task->SetFillBayesianPIDInfo(kFALSE);
  
  // Event and track selection
  //--------------------------
  task->SetEventFilter(CreateEventFilter(isAOD));
  task->SetTrackFilter(CreateGlobalTrackFilter(isAOD));
  //  task->SetFlowTrackFilter(CreateFlowTrackFilter(isAOD));
  //task->SetK0sPionCuts(CreateK0sPionCuts(isAOD));
  //task->SetLambdaProtonCuts(CreateLambdaProtonCuts(isAOD));
  //task->SetLambdaPionCuts(CreateLambdaPionCuts(isAOD));
  //task->SetGammaElectronCuts(CreateGammaConvElectronCuts(isAOD));
  //task->SetK0sMassRange(0.44,0.55);
  //task->SetLambdaMassRange(1.090,1.14);
  //task->SetGammaConvMassRange(0.0,0.1);
  //task->SetV0OpenCuts(CreateV0OpenCuts(AliESDv0KineCuts::kPurity, AliESDv0KineCuts::kPbPb));
  //task->SetV0StrongCuts(CreateV0StrongCuts(AliESDv0KineCuts::kPurity, AliESDv0KineCuts::kPbPb));
  task->SetFillFMDInfo(kFALSE);
  //AddFMDTask();
  //NEW: for Reaction plane:
  task->SetFillEventPlaneInfo(kFALSE);

  // Set Writting options
  //---------------------
  //task->SetTreeWritingOption(AliAnalysisTaskReducedTreeMaker::kFullEventsWithFullTracks);
  //task->SetTreeWritingOption(AliAnalysisTaskReducedTreeMaker::kBaseEventsWithBaseTracks);
  if(reducedEventType!=-1)
    task->SetTreeWritingOption(reducedEventType);
  task->SetWriteTree(writeTree);

  // For MC:
  //--------
  task->SetFillMCInfo(kTRUE);

  // You can define the active branches of the tree, all other branches will be set to inactive
  //-------------------------------------------------------------------------------------------
  //task->SetTreeActiveBranch("Event");
  //task->SetTreeActiveBranch("fRunNo");
  // Alternatively, you can define the inactive branches of the tree, all other branches will be set to active
  task->SetTreeInactiveBranch("fUniqueID");
  task->SetTreeInactiveBranch("fBits");
  task->SetTreeInactiveBranch("fEventTag");
  task->SetTreeInactiveBranch("fNVtxContributors");
  task->SetTreeInactiveBranch("fCentQuality");
  task->SetTreeInactiveBranch("fNtracks*");
  task->SetTreeInactiveBranch("fNV0candidates*");
  task->SetTreeInactiveBranch("fTracks.fUniqueID");
  task->SetTreeInactiveBranch("fTracks.fBits");
  task->SetTreeInactiveBranch("fTracks.fIsCartesian");
  task->SetTreeInactiveBranch("fTracks.fFlags");
  task->SetTreeInactiveBranch("fTracks.fQualityFlags");
  task->SetTreeInactiveBranch("fTracks.fTrackId");
  task->SetTreeInactiveBranch("fTracks.fStatus");
  task->SetTreeInactiveBranch("fTracks.fTPCPhi");
  task->SetTreeInactiveBranch("fTracks.fTPCPt");
  task->SetTreeInactiveBranch("fTracks.fTPCEta");
  task->SetTreeInactiveBranch("fTracks.fMomentumInner");
  //task->SetTreeInactiveBranch("fTracks.fDCA*");
  task->SetTreeInactiveBranch("fTracks.fTPCDCA*");
  task->SetTreeInactiveBranch("fTracks.fTrackLength");
  task->SetTreeInactiveBranch("fTracks.fMassForTracking");
  task->SetTreeInactiveBranch("fTracks.fChi2TPCConstrainedVsGlobal");
  //task->SetTreeInactiveBranch("fTracks.fITSclusterMap");
  task->SetTreeInactiveBranch("fTracks.fITSSharedclusterMap");
  task->SetTreeInactiveBranch("fTracks.fITSsignal");
  task->SetTreeInactiveBranch("fTracks.fITSchi2");
  task->SetTreeInactiveBranch("fTracks.fTPCNcls");
  task->SetTreeInactiveBranch("fTracks.fTPCCrossedRows");
  task->SetTreeInactiveBranch("fTracks.fTPCNclsF");
  task->SetTreeInactiveBranch("fTracks.fTPCNclsShared");
  task->SetTreeInactiveBranch("fTracks.fTPCClusterMap");
  task->SetTreeInactiveBranch("fTracks.fTPCsignal");
  task->SetTreeInactiveBranch("fTracks.fTPCsignalN");
  task->SetTreeInactiveBranch("fTracks.fTPCnSig*");
  task->SetTreeInactiveBranch("fTracks.fTPCchi2");
  task->SetTreeInactiveBranch("fTracks.fTPCActiveLength");
  task->SetTreeInactiveBranch("fTracks.fTPCGeomLength");
  task->SetTreeInactiveBranch("fTracks.fTOFbeta");
  task->SetTreeInactiveBranch("fTracks.fTOFmismatchProbab");
  task->SetTreeInactiveBranch("fTracks.fTRDntracklets*");
  task->SetTreeInactiveBranch("fTracks.fMCMom*");
  task->SetTreeInactiveBranch("fTracks.fMCFreezeout*");
  task->SetTreeInactiveBranch("fTracks.fMCLabels*");
  task->SetTreeInactiveBranch("fTracks.fMCPdg*");
  task->SetTreeInactiveBranch("fTracks.fMCGeneratorIndex");
  task->SetTreeInactiveBranch("fTracks.fTOFtime");
  task->SetTreeInactiveBranch("fTracks.fTOFdx");
  task->SetTreeInactiveBranch("fTracks.fTOFdz");
  task->SetTreeInactiveBranch("fTracks.fTOFchi2");
  task->SetTreeInactiveBranch("fTracks.fTOFnSig*");
  task->SetTreeInactiveBranch("fTracks.fTOFdeltaBC");
  task->SetTreeInactiveBranch("fTracks.fTRDpid*");
  task->SetTreeInactiveBranch("fTracks.fTRDpidLQ2D*");
  task->SetTreeInactiveBranch("fTracks.fCaloClusterId");
  task->SetTreeInactiveBranch("fCandidates.*");
  task->SetTreeInactiveBranch("fEventNumberInFile");
  //task->SetTreeInactiveBranch("fIsPhysicsSelection");
  task->SetTreeInactiveBranch("fVtxSPD*");
  task->SetTreeInactiveBranch("fNVtxSPDContributors");
  task->SetTreeInactiveBranch("fIsSPDPileup");
  task->SetTreeInactiveBranch("fIsSPDPileupMultBins");
  task->SetTreeInactiveBranch("fSPDntracklets");
  task->SetTreeInactiveBranch("fITSClusters*");
  task->SetTreeInactiveBranch("fSPDnSingle");
  task->SetTreeInactiveBranch("fL0TriggerInputs");
  task->SetTreeInactiveBranch("fL1TriggerInputs");
  task->SetTreeInactiveBranch("fL2TriggerInputs");
  task->SetTreeInactiveBranch("fBC");
  task->SetTreeInactiveBranch("fTimeStamp");
  task->SetTreeInactiveBranch("fEventType");
  task->SetTreeInactiveBranch("fTriggerMask");
  task->SetTreeInactiveBranch("fIRIntClosestIntMap*");
  task->SetTreeInactiveBranch("fVtxTPC*");
  task->SetTreeInactiveBranch("fNVtxTPCContributors");
  task->SetTreeInactiveBranch("fNpileupSPD");
  task->SetTreeInactiveBranch("fNpileupTracks");
  task->SetTreeInactiveBranch("fNPMDtracks");
  task->SetTreeInactiveBranch("fNTRDtracks");
  task->SetTreeInactiveBranch("fNTRDtracklets");
  task->SetTreeInactiveBranch("fSPDntrackletsEta*");
  task->SetTreeInactiveBranch("fSPDFiredChips*");
  task->SetTreeInactiveBranch("fNtracksPerTrackingFlag*");
  task->SetTreeInactiveBranch("fVZEROMult*");
  task->SetTreeInactiveBranch("fVZEROTotalMult*");
  task->SetTreeInactiveBranch("fZDCnTotalEnergy*");
  task->SetTreeInactiveBranch("fZDCpTotalEnergy*");
  task->SetTreeInactiveBranch("fZDCnEnergy*");
  task->SetTreeInactiveBranch("fZDCpEnergy*");
  task->SetTreeInactiveBranch("fT0amplitude*");
  task->SetTreeInactiveBranch("fT0TOF*");
  task->SetTreeInactiveBranch("fT0TOFbest*");
  task->SetTreeInactiveBranch("fT0zVertex");
  task->SetTreeInactiveBranch("fT0start");
  task->SetTreeInactiveBranch("fT0pileup");
  task->SetTreeInactiveBranch("fT0sattelite");
  task->SetTreeInactiveBranch("fNCaloClusters");
  
  // Add task to the manager
  //------------------------
  mgr->AddTask(task);
    
  // Output containers
  //------------------
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

  // General track cuts
  //-------------------

  // -- cuts group
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
 
  // -- Recommended Pb-Pb runII cuts:
  //static AliESDtrackCuts* GetStandardITSTPCTrackCuts2015PbPb(Bool_t selPrimaries=kTRUE, Int_t clusterCut=1, Bool_t cutAcceptanceEdges = kTRUE, Bool_t removeDistortedRegions = kFALSE);
  AliESDtrackCuts* esdTrackCuts = 0x0;
  esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb(kTRUE, 1, kTRUE, kFALSE);

  // -- Additional cuts:
  AliESDtrackCuts *esdTrackCuts2 = new AliESDtrackCuts;
  esdTrackCuts2->SetPtRange( 0.2 , 200.0);
  //esdTrackCuts2->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
  esdTrackCuts2->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  cuts->AddCut(esdTrackCuts2);

  cuts->AddCut(esdTrackCuts);

  // PID cuts
  //---------

  //pass2_lowIR for production (TOF-IF) pT>0.2GeV:
  AliDielectronPID *ePID = new AliDielectronPID("noTOF","noTOF");
  ePID->AddCut(AliDielectronPID::kITS, AliPID::kElectron,     -3., 3.,     0.2, 1.5,    kFALSE );
  ePID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron ,    -2., 2.,     0.4, 100.,   kFALSE );
  ePID->AddCut(AliDielectronPID::kTPC, AliPID::kElectron ,    -3., 3.,     0.2,  0.4,   kFALSE );
  ePID->AddCut(AliDielectronPID::kTPC, AliPID::kPion,       -100., 4.,     0.2, 100.,   kTRUE  );
  ePID->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,     -3., 3.,     0.2,   5.,   kFALSE , AliDielectronPID::kIfAvailable );

  cuts->AddCut(ePID);
 
  return cuts;
}
