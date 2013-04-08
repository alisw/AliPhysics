AliAnalysisTaskMinijet* AddTaskMinijet(Int_t runNumber     =    -1, 
				       TString format      = "aod",
				       Bool_t  useMC       = false, 
				       Bool_t  mcOnly      = false,
				       Bool_t  useHighMult = false,
				       Float_t ptTrigMin   =   0.7,
				       Float_t ptAssocMin  =   0.4,
				       Float_t maxVtxZ     =   10.,
				       Int_t   filterBit   =    96,
				       Int_t   debugLevel  =     0,
				       Float_t maxEta      =   0.9,
				       Float_t ptMin       =   0.2,
				       Float_t ptMax       =  50.0,
				       Bool_t  checkSDD    = false,
				       Int_t   checkOption =     1,
				       Bool_t  onlyPrim    = true,
				       Bool_t  corrStrange = true,
				       Bool_t  threePart   = false,
				       Bool_t  rejectChunk = false,
				       Int_t   nTPC        = 5,
				       const char *foldername = "MiniJets")
{
  
  // Get the pointer to the existing analysis manager via the static access method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("Analysis train");

  // Set TPC track cuts (used for ESDs) 
  AliESDtrackCuts* esdTrackCutsTPC=0x0;
  if(!format.CompareTo("esd")){
    esdTrackCutsTPC = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    esdTrackCutsTPC->SetMinNClustersTPC(80);
  }

  // First task for min bias events
  AliAnalysisTaskMinijet *taskMB = 
    new AliAnalysisTaskMinijet("AliAnalysisTaskMinijet Min bias");
  taskMB->SetUseMC(useMC,mcOnly);
  taskMB->SetTriggerPtCut(ptTrigMin);
  taskMB->SetAssociatePtCut(ptAssocMin);
  taskMB->SetMaxEta(maxEta);
  taskMB->SetPtRange(ptMin, ptMax);
  taskMB->SetMaxVertexZ(maxVtxZ);
  taskMB->SetTriggerMask(AliVEvent::kMB);
  // taskMB->SelectCollisionCandidates(AliVEvent::kMB);//MB  //now inside task
  taskMB->SetFilterBit(filterBit); // used only in AOD case
  taskMB->SetCheckSDD(checkSDD, checkOption);
  taskMB->SetCorrStrangeness(corrStrange);
  taskMB->SetThreeParticleCorrelation(threePart);
  taskMB->SetRejectCorrupted(rejectChunk, nTPC);
  taskMB->SetDebugLevel(debugLevel);

  //use this only for correction map plots -> Split contmaintion and rec efficiency
  taskMB->SetAnalyseOnlyPrimaries(onlyPrim);
 
  if(!format.CompareTo("esd")){
    taskMB->SetCuts(esdTrackCutsTPC);
    taskMB->SetModeEsdAod(0); // 0 = reading ESDs
  }  
  else if (!format.CompareTo("aod")){
    // Cuts already applied through esd filter
    taskMB->SetModeEsdAod(1); // 1 = reading AODs
  }

  // Second task for high multipliciy events
  AliAnalysisTaskMinijet *taskHM =0x0;
  if(useHighMult){
    taskHM  = new AliAnalysisTaskMinijet("AliAnalysisTaskMinijet HighMult");
    taskHM->SetUseMC(useMC, mcOnly);
    taskHM->SetTriggerPtCut(ptTrigMin);
    taskHM->SetAssociatePtCut(ptAssocMin);
    taskHM->SetMaxEta(maxEta);
    taskHM->SetPtRange(ptMin, ptMax);
    taskHM->SetMaxVertexZ(maxVtxZ);
    taskHM->SetTriggerMask(AliVEvent::kHighMult);
    //taskHM->SelectCollisionCandidates(AliVEvent::kHighMult); // now inside task
    taskHM->SetFilterBit(filterBit); // used only in AOD case
    taskHM->SetCheckSDD(checkSDD, checkOption);
    taskHM->SetCorrStrangeness(corrStrange);
    taskHM->SetThreeParticleCorrelation(threePart);
    taskHM->SetRejectCorrupted(rejectChunk, nTPC);
    taskHM->SetDebugLevel(debugLevel);

    if(!format.CompareTo("esd")){
      taskHM->SetCuts(esdTrackCutsTPC);
      taskHM->SetModeEsdAod(0); // 0 = reading ESDs
    }
    else if (!format.CompareTo("aod")){
      //cuts already applied through esd filter during writing of AODs
      taskHM->SetModeEsdAod(1); // 1 = reading AODs
    }
    
  }

  //create output containers 
  AliAnalysisDataContainer *outputMB  = 0x0;
  AliAnalysisDataContainer *outputHM = 0x0;

  if(runNumber>0){ 
    outputMB  =  mgr->CreateContainer(foldername,TList::Class(),
				      AliAnalysisManager::kOutputContainer, 
				      Form("run%d.root",runNumber));
    if(useHighMult){
      outputHM =  mgr->CreateContainer("MiniJets_HighMult",TList::Class(),
				       AliAnalysisManager::kOutputContainer, 
				       Form("run%d.root",runNumber));
    }
  }
  else{
    outputMB  = mgr->CreateContainer(foldername,TList::Class(),
				     AliAnalysisManager::kOutputContainer, 
				     Form("%s:PWG4_MiniJets",
					  AliAnalysisManager::GetCommonFileName()));
    if(useHighMult){
      outputHM = mgr->CreateContainer("MiniJets_HighMult",TList::Class(),
				      AliAnalysisManager::kOutputContainer, 
				      Form("%s:PWG4_MiniJets",
					   AliAnalysisManager::GetCommonFileName()));
    }
  }
  

  // Add first task to the manager and connect container
  mgr->AddTask(taskMB);
  mgr->ConnectInput(taskMB, 0,  mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskMB, 1, outputMB);
  
  
  // Add second task to the manager and connect container
  if(useHighMult){
    mgr->AddTask(taskHM);
    mgr->ConnectInput (taskHM, 0,  mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskHM, 1, outputHM);
  }
  
  return taskMB;

}

