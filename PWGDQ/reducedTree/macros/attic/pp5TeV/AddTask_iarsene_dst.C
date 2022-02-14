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
  if(!prod.CompareTo("LHC15n") || !prod.CompareTo("LHC17p") || !prod.CompareTo("LHC17q")) collSystem = 9;
  if(!prod.CompareTo("LHC17k") || !prod.CompareTo("LHC17o") || !prod.CompareTo("LHC18c")) collSystem = 10;
  
  if(!collSystem && !hasMC) {
     collSystem = 1;
     printf("WARNING: In AddTask_iarsene_dst(), no proper production name specified, or not supported! \n");
     printf("                 Using collSystem=1 (min bias 2010 Pb-Pb) as default \n");
  }
  
  //create task and add it to the manager
  AliAnalysisTaskReducedTreeMaker *task=new AliAnalysisTaskReducedTreeMaker("DSTTreeMaker", kTRUE);
  if(collSystem==1) task->SetTriggerMask(AliVEvent::kMB);
  if(collSystem==2)
    task->SetTriggerMask(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  if(collSystem==3) task->SetTriggerMask(AliVEvent::kINT7);
  if(collSystem==4) task->SetTriggerMask(AliVEvent::kMB);
  if(collSystem==5) task->SetTriggerMask(AliVEvent::kINT7);
  if(collSystem==6) task->SetTriggerMask(AliVEvent::kMB | AliVEvent::kINT7 | AliVEvent::kHighMultSPD | AliVEvent::kHighMultV0);
  if(collSystem==7) task->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kTRD);
  if(collSystem==10) task->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kTRD | AliVEvent::kHighMultV0 | AliVEvent::kEMCEJE | AliVEvent::kEMCEGA);
  //if(collSystem==10) task->SetTriggerMask(AliVEvent::kINT7);
  if(collSystem==8) task->SetTriggerMask(AliVEvent::kINT7);
  if(collSystem==9) task->SetTriggerMask(AliVEvent::kINT7);
  if(collSystem==0) task->SetTriggerMask(AliVEvent::kINT7);
  
  //if(collSystem==4)
  //task->SetRejectPileup();
    
  task->UsePhysicsSelection(kTRUE);
  task->SetUseAnalysisUtils(kTRUE);
  
  task->SetFillV0Info(kFALSE);
  task->SetFillGammaConversions(kFALSE);
  task->SetFillK0s(kFALSE);
  task->SetFillLambda(kFALSE);
  task->SetFillALambda(kFALSE);
  task->SetFillCaloClusterInfo(kFALSE);
  //task->SetFillDielectronInfo(kFALSE);
  //task->SetFillFriendInfo(kFALSE);
  
  task->SetEventFilter(CreateEventFilter(isAOD));
  //task->SetTrackFilter(CreateGlobalTrackFilter(isAOD));
  CreateTrackFilters(task, isAOD);
  
  //  task->SetFlowTrackFilter(CreateFlowTrackFilter(isAOD));
  //task->SetK0sPionCuts(CreateK0sPionCuts(isAOD));
  //task->SetLambdaProtonCuts(CreateLambdaProtonCuts(isAOD));
  //task->SetLambdaPionCuts(CreateLambdaPionCuts(isAOD));
  task->SetGammaElectronCuts(CreateGammaConvElectronCuts(isAOD));
  //task->SetK0sMassRange(0.44,0.55);
  //task->SetLambdaMassRange(1.090,1.14);
  task->SetGammaConvMassRange(0.0,0.1);
  task->SetV0OpenCuts(CreateV0OpenCuts(AliESDv0KineCuts::kPurity, AliESDv0KineCuts::kPbPb));
  task->SetV0StrongCuts(CreateV0StrongCuts(AliESDv0KineCuts::kPurity, AliESDv0KineCuts::kPbPb));
  //task->SetFillFMDInfo(); AddFMDTask();
  if(hasMC) {
     task->SetFillMCInfo(kTRUE);
     AddMCSignals(task);
  }
  //task->SetFillEventPlaneInfo(kTRUE);

  //task->SetTreeWritingOption(AliAnalysisTaskReducedTreeMaker::kFullEventsWithFullTracks);
  //task->SetTreeWritingOption(AliAnalysisTaskReducedTreeMaker::kBaseEventsWithBaseTracks);
  if(reducedEventType!=-1)
    task->SetTreeWritingOption(reducedEventType);
  task->SetWriteTree(writeTree);
  //if(!hasMC)
   //task->SetEventWritingRequirement(2,0,0.005);
  //task->SetWriteEventsWithNoSelectedTracks(kFALSE, 0.05, 2);
  //task->SetWriteEventsWithNoSelectedTracksAndNoSelectedAssociatedTracks(kFALSE);
  
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
   
   AliSignalMC* jpsiFromB=new AliSignalMC("JpsiNonPrompt","",1,2);
   jpsiFromB->SetPDGcode(0, 0, 443, kFALSE);
   jpsiFromB->SetPDGcode(0, 1, 503, kTRUE);
   task->AddMCsignal(jpsiFromB, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   
   AliSignalMC* jpsiPrompt=new AliSignalMC("JpsiPrompt","",1,2);
   jpsiPrompt->SetPDGcode(0, 0, 443, kFALSE);
   jpsiPrompt->SetPDGcode(0, 1, 503, kTRUE, kTRUE);
   task->AddMCsignal(jpsiPrompt, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   
   AliSignalMC* jpsiInclusiveRadiative=new AliSignalMC("JpsiInclusiveRadiative","",1,1);
   jpsiInclusiveRadiative->SetPDGcode(0, 0, 443, kFALSE);
   jpsiInclusiveRadiative->SetSourceBit(0, 0, AliSignalMC::kRadiativeDecay);
   task->AddMCsignal(jpsiInclusiveRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   
   AliSignalMC* jpsiInclusiveNonRadiative=new AliSignalMC("JpsiInclusiveNonRadiative","",1,1);
   jpsiInclusiveNonRadiative->SetPDGcode(0, 0, 443, kFALSE);
   jpsiInclusiveNonRadiative->SetSourceBit(0, 0, AliSignalMC::kRadiativeDecay, kTRUE);
   task->AddMCsignal(jpsiInclusiveNonRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   
   AliSignalMC* electronFromJpsi=new AliSignalMC("electronFromJpsiInclusive","",1,2);
   electronFromJpsi->SetPDGcode(0, 0, 11, kTRUE);
   electronFromJpsi->SetPDGcode(0, 1, 443);
   task->AddMCsignal(electronFromJpsi, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   
   AliSignalMC* electronFromJpsiNonPrompt=new AliSignalMC("electronFromJpsiNonPrompt","",1,3);
   electronFromJpsiNonPrompt->SetPDGcode(0, 0, 11, kTRUE);
   electronFromJpsiNonPrompt->SetPDGcode(0, 1, 443);
   electronFromJpsiNonPrompt->SetPDGcode(0, 2, 503, kTRUE);      //
   task->AddMCsignal(electronFromJpsiNonPrompt, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   
   AliSignalMC* electronFromJpsiPrompt=new AliSignalMC("electronFromJpsiPrompt","", 1,3);
   electronFromJpsiPrompt->SetPDGcode(0, 0, 11, kTRUE);
   electronFromJpsiPrompt->SetPDGcode(0, 1, 443);
   electronFromJpsiPrompt->SetPDGcode(0, 2, 503, kTRUE, kTRUE);    //
   task->AddMCsignal(electronFromJpsiPrompt, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   
   AliSignalMC* electronFromJpsiRadiative=new AliSignalMC("electronFromJpsiRadiative","",1,2);
   electronFromJpsiRadiative->SetPDGcode(0, 0, 11, kTRUE);
   electronFromJpsiRadiative->SetPDGcode(0, 1, 443);
   electronFromJpsiRadiative->SetSourceBit(0, 1, AliSignalMC::kRadiativeDecay);
   task->AddMCsignal(electronFromJpsiRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   
   AliSignalMC* electronFromJpsiNonRadiative=new AliSignalMC("electronFromJpsiNonRadiative","",1,2);
   electronFromJpsiNonRadiative->SetPDGcode(0, 0, 11, kTRUE);
   electronFromJpsiNonRadiative->SetPDGcode(0, 1, 443);
   electronFromJpsiNonRadiative->SetSourceBit(0, 1, AliSignalMC::kRadiativeDecay, kTRUE);
   task->AddMCsignal(electronFromJpsiNonRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   
   AliSignalMC* photonFromJpsiDecay=new AliSignalMC("photonFromJpsiDecay","",1,2);
   photonFromJpsiDecay->SetPDGcode(0, 0, 22);
   photonFromJpsiDecay->SetPDGcode(0, 1, 443);
   task->AddMCsignal(photonFromJpsiDecay, AliAnalysisTaskReducedTreeMaker::kFullTrack);
   
   AliSignalMC* primary = new AliSignalMC("primary", "", 1, 1);
   primary->SetSourceBit(0, 0, AliSignalMC::kPhysicalPrimary);
   //task->AddMCsignal(primary, AliAnalysisTaskReducedTreeMaker::kFullTrack);
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
   //task->SetTreeInactiveBranch("fNV0candidates*");
   //task->SetTreeInactiveBranch("fTracks.*");
   task->SetTreeInactiveBranch("fCandidates.*");
   task->SetTreeInactiveBranch("fEventNumberInFile");
   task->SetTreeInactiveBranch("fL0TriggerInputs");
   task->SetTreeInactiveBranch("fL1TriggerInputs");
   task->SetTreeInactiveBranch("fL2TriggerInputs");
   //task->SetTreeInactiveBranch("fBC");
   //task->SetTreeInactiveBranch("fTimeStamp");
   task->SetTreeInactiveBranch("fEventType");
   //task->SetTreeInactiveBranch("fTriggerMask");
   //task->SetTreeInactiveBranch("fMultiplicityEstimators*");
   //task->SetTreeInactiveBranch("fMultiplicityEstimatorPercentiles*");
   //task->SetTreeInactiveBranch("fIsPhysicsSelection");
   //task->SetTreeInactiveBranch("fIsSPDPileup");
   //task->SetTreeInactiveBranch("fIsSPDPileupMultBins");
   //task->SetTreeInactiveBranch("fIRIntClosestIntMap*");
   //task->SetTreeInactiveBranch("fVtxTPC*");
   task->SetTreeInactiveBranch("fVtxMC*");
   task->SetTreeInactiveBranch("fNch*");
   //task->SetTreeInactiveBranch("fNVtxTPCContributors");
   //task->SetTreeInactiveBranch("fVtxSPD*");
   //task->SetTreeInactiveBranch("fNVtxSPDContributors");
   //task->SetTreeInactiveBranch("fNpileupSPD");
   //task->SetTreeInactiveBranch("fNpileupTracks");
   //task->SetTreeInactiveBranch("fNTPCclusters");
   task->SetTreeInactiveBranch("fNPMDtracks");
   task->SetTreeInactiveBranch("fNTRDtracks");
   task->SetTreeInactiveBranch("fNTRDtracklets");
   //task->SetTreeInactiveBranch("fSPDntracklets");
   //task->SetTreeInactiveBranch("fSPDntrackletsEta*");
   //task->SetTreeInactiveBranch("fSPDFiredChips*");
   //task->SetTreeInactiveBranch("fITSClusters*");
   //task->SetTreeInactiveBranch("fSPDnSingle");
   //task->SetTreeInactiveBranch("fNtracksPerTrackingFlag*");
   //task->SetTreeInactiveBranch("fVZEROMult*");
   //task->SetTreeInactiveBranch("fVZEROTotalMult*");
   task->SetTreeInactiveBranch("fZDCnEnergy*");
   task->SetTreeInactiveBranch("fZDCpEnergy*");
   //task->SetTreeInactiveBranch("fZDCnTotalEnergy*");
   //task->SetTreeInactiveBranch("fZDCpTotalEnergy*");
   task->SetTreeInactiveBranch("fT0amplitude*");
   task->SetTreeInactiveBranch("fT0TOF*");
   task->SetTreeInactiveBranch("fT0TOFbest*");
   task->SetTreeInactiveBranch("fT0zVertex");
   task->SetTreeInactiveBranch("fT0start");
   //task->SetTreeInactiveBranch("fT0pileup");
   //task->SetTreeInactiveBranch("fT0sattelite");
   //task->SetTreeInactiveBranch("fNCaloClusters");
   //task->SetTreeInactiveBranch("fCaloClusters.*");
   //task->SetTreeInactiveBranch("fFMD.*");
   //task->SetTreeInactiveBranch("fEventPlane.*");
   
   //task->SetTreeInactiveBranch("fTracks2.*");
   //task->SetTreeInactiveBranch("fTracks.fP*");
   //task->SetTreeInactiveBranch("fTracks.fIsCartesian");
   //task->SetTreeInactiveBranch("fTracks.fCharge");
   //task->SetTreeInactiveBranch("fTracks.fFlags");
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
   //task->SetTreeInactiveBranch("fTracks.fTOFmismatchProbab");
   //task->SetTreeInactiveBranch("fTracks.fTOFchi2");
   //task->SetTreeInactiveBranch("fTracks.fTOFnSig*");
   //task->SetTreeInactiveBranch("fTracks.fTOFdeltaBC");
   task->SetTreeInactiveBranch("fTracks.fTRDntracklets*");
   task->SetTreeInactiveBranch("fTracks.fTRDpid*");
   task->SetTreeInactiveBranch("fTracks.fTRDpidLQ2D*");
   task->SetTreeInactiveBranch("fTracks.fCaloClusterId");
   //task->SetTreeInactiveBranch("fTracks.fMCMom*");
   //task->SetTreeInactiveBranch("fTracks.fMCFreezeout*");
   //task->SetTreeInactiveBranch("fTracks.fMCLabels*");
   //task->SetTreeInactiveBranch("fTracks.fMCPdg*");
   //task->SetTreeInactiveBranch("fTracks.fMCGeneratorIndex");
  
   task->SetTreeInactiveBranch("fCandidates.*");
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
  //if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  //eventCuts->SetRequireVertex();
  //eventCuts->SetMinVtxContributors(1);
  //eventCuts->SetVertexZ(-10.,10.);
  return eventCuts;
}


//______________________________________________________________________________________
void CreateTrackFilters(AliAnalysisTaskReducedTreeMaker* task, Bool_t isAOD) {
   //
   // add track filters
   //
   AliDielectronCutGroup* basicCut = new AliDielectronCutGroup("basicCut","J/psi candidate electron (basic cut)",AliDielectronCutGroup::kCompAND);
   AliDielectronVarCuts *basicVarCuts = new AliDielectronVarCuts("basicVarCuts","basic var track cuts");
   basicVarCuts->AddCut(AliDielectronVarManager::kImpactParXY,-1.0,1.0);
   basicVarCuts->AddCut(AliDielectronVarManager::kImpactParZ,-3.0,3.0);
   basicVarCuts->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
   basicVarCuts->AddCut(AliDielectronVarManager::kP,0.9,1.0e+30);
   basicVarCuts->AddCut(AliDielectronVarManager::kNclsTPC,70.0,161.0);   
   basicVarCuts->AddCut(AliDielectronVarManager::kTPCnSigmaEle, -3.0, 4.0);
   basicVarCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPro, -3.0, 2.0, kTRUE);
   //basicVarCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPio, -3.0, 2.0, kTRUE);
   basicVarCuts->AddCut(AliDielectronVarManager::kTPCnSigmaKao, -3.0, 3.0, kTRUE);
//   basicVarCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl, 0.0, 10.0);
//   basicVarCuts->AddCut(AliDielectronVarManager::kITSchi2Cl, 0.0, 100.);
//   basicVarCuts->AddCut(AliDielectronVarManager::kTPCclsSegments, 5., 9.);
//   basicVarCuts->AddCut(AliDielectronVarManager::kNclsSFracTPC, 0.6, 2.0, kTRUE);
//   basicVarCuts->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.5, 2.0);
//   basicVarCuts->AddCut(AliDielectronVarManager::kChi2TPCConstrainedVsGlobal, 100., 1.0e+10, kTRUE);
   basicCut->AddCut(basicVarCuts);
   AliDielectronTrackCuts* basicTrackCuts = new AliDielectronTrackCuts("basicTrackCuts","basic track cuts");
   basicTrackCuts->SetRequireITSRefit(kTRUE);
   basicTrackCuts->SetRequireTPCRefit(kTRUE);
   if(isAOD) basicTrackCuts->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
   basicCut->AddCut(basicTrackCuts);
   task->AddTrackFilter(basicCut, kFALSE);   // basic cut, no SPD requirements
   
   AliDielectronCutGroup* prefilterCut = new AliDielectronCutGroup("prefilterCut","basic electron cut to be used for prefilter",AliDielectronCutGroup::kCompAND);
   AliDielectronVarCuts *prefilterVarCuts = new AliDielectronVarCuts("prefilterVarCuts","prefilter var track cuts");
   prefilterVarCuts->AddCut(AliDielectronVarManager::kImpactParXY,-3.0,3.0);
   prefilterVarCuts->AddCut(AliDielectronVarManager::kImpactParZ,-10.0,10.0);
   prefilterVarCuts->AddCut(AliDielectronVarManager::kEta,-0.9,0.9);
   prefilterVarCuts->AddCut(AliDielectronVarManager::kPt,0.2,1.0e+30);
   prefilterVarCuts->AddCut(AliDielectronVarManager::kNclsTPC,70.0,161.0);   
   prefilterVarCuts->AddCut(AliDielectronVarManager::kTPCnSigmaEle, -4.0, 4.0);
   prefilterVarCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPro, -3.0, 3.0, kTRUE);
   //prefilterVarCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPio, -3.0, 3.0, kTRUE);
   prefilterVarCuts->AddCut(AliDielectronVarManager::kTPCnSigmaKao, -3.0, 3.0, kTRUE);
//   prefilterVarCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl, 0.0, 10.0);
//   prefilterVarCuts->AddCut(AliDielectronVarManager::kITSchi2Cl, 0.0, 100.);
//   prefilterVarCuts->AddCut(AliDielectronVarManager::kTPCclsSegments, 4., 9.);
   //prefilterVarCuts->AddCut(AliDielectronVarManager::kNclsSFracTPC, 0.6, 2.0, kTRUE);
   //prefilterVarCuts->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.5, 2.0);
   //prefilterVarCuts->AddCut(AliDielectronVarManager::kChi2TPCConstrainedVsGlobal, 100., 1.0e+10, kTRUE);
   prefilterCut->AddCut(prefilterVarCuts);
   AliDielectronTrackCuts* prefilterTrackCuts = new AliDielectronTrackCuts("prefilterTrackCuts","prefilter track cuts");
   //prefilterTrackCuts->SetRequireITSRefit(kTRUE);
   prefilterTrackCuts->SetRequireTPCRefit(kTRUE);
   if(isAOD) prefilterTrackCuts->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
   prefilterCut->AddCut(prefilterTrackCuts);
   task->AddTrackFilter(prefilterCut, kTRUE);   // prefilter electron cut
   
   AliDielectronCutGroup* assocHadr = new AliDielectronCutGroup("assocHadr","Hadrons",AliDielectronCutGroup::kCompAND);
   AliDielectronVarCuts *hadrVarCuts = new AliDielectronVarCuts("hadrVarCuts","hadron var track cuts");
   hadrVarCuts->AddCut(AliDielectronVarManager::kImpactParXY,-0.5,0.5);
   hadrVarCuts->AddCut(AliDielectronVarManager::kImpactParZ,-1.0,1.0);
   hadrVarCuts->AddCut(AliDielectronVarManager::kEta,-1.0,1.0);
   hadrVarCuts->AddCut(AliDielectronVarManager::kPt,0.2,1.0e+30);
   hadrVarCuts->AddCut(AliDielectronVarManager::kNclsTPC,70.0,161.0);   
   hadrVarCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl, 0.0, 4.0);
   hadrVarCuts->AddCut(AliDielectronVarManager::kITSchi2Cl, 0.0, 20.);
   assocHadr->AddCut(hadrVarCuts);
   assocHadr->AddCut(basicTrackCuts);     // ITSrefit and TPCrefit
   //task->AddTrackFilter(assocHadr, kTRUE);
   
   AliDielectronCutGroup* assocHadr2 = new AliDielectronCutGroup("assocHadrNoITSrefit","Hadrons, no ITSrefit",AliDielectronCutGroup::kCompAND);
   assocHadr2->AddCut(hadrVarCuts);
   assocHadr2->AddCut(prefilterTrackCuts);    // TPCrefit (no ITSrefit  --> check consistency !!!!!!!!!!!!)
   //task->AddTrackFilter(assocHadr2, kTRUE);
   
   AliDielectronCutGroup* assocPion = new AliDielectronCutGroup("assocPion","Associated pions",AliDielectronCutGroup::kCompAND);
   assocPion->AddCut(hadrVarCuts);
   assocPion->AddCut(basicTrackCuts);
   AliDielectronVarCuts *pionPID = new AliDielectronVarCuts("pionPID","pion PID");
   pionPID->AddCut(AliDielectronVarManager::kTPCnSigmaPio, -3.0, 3.0);
   assocPion->AddCut(pionPID);
   //task->AddTrackFilter(assocPion, kTRUE);
   
   AliDielectronCutGroup* assocKaon = new AliDielectronCutGroup("assocKaon","Associated kaons",AliDielectronCutGroup::kCompAND);
   assocKaon->AddCut(hadrVarCuts);
   assocKaon->AddCut(basicTrackCuts);
   AliDielectronVarCuts *kaonPID = new AliDielectronVarCuts("kaonPID","kaon PID");
   kaonPID->AddCut(AliDielectronVarManager::kTPCnSigmaKao, -3.0, 3.0);
   assocKaon->AddCut(kaonPID);
   //task->AddTrackFilter(assocKaon, kTRUE);
   
   AliDielectronCutGroup* assocProton = new AliDielectronCutGroup("assocProton","Associated protons",AliDielectronCutGroup::kCompAND);
   assocProton->AddCut(hadrVarCuts);
   assocProton->AddCut(basicTrackCuts);
   AliDielectronVarCuts *protonPID = new AliDielectronVarCuts("protonPID","proton PID");
   protonPID->AddCut(AliDielectronVarManager::kTPCnSigmaPro, -3.0, 3.0);
   assocProton->AddCut(protonPID);
   //task->AddTrackFilter(assocProton, kTRUE);
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
    electronCuts->SetPtRange(0.6,100.0);
    electronCuts->SetEtaRange(-0.9,0.9);
    electronCuts->SetRequireTPCRefit(kTRUE);
    electronCuts->SetMinNClustersTPC(70);
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
