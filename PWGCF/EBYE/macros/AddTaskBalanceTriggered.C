// now in options

//PID config
Bool_t kUseNSigmaPID = kFALSE;
Double_t nSigmaMax = 3.0;
Bool_t kUseBayesianPID = kTRUE;
Double_t gMinAcceptedProbability = 0.7;

//_________________________________________________________//
AliAnalysisTaskTriggeredBF *AddTaskBalanceTriggered(Double_t centrMin=0.,
						 Double_t centrMax=100.,
						 Bool_t gRunShuffling=kFALSE,
						 Bool_t gRunMixing=kFALSE,
						 TString centralityEstimator="V0M",
						 Double_t vertexZ=10.,
						 Double_t DCAxy=-1,
						 Double_t DCAz=-1,
						 Double_t ptMin=0.15,
						 Double_t ptMax=20,
						 Double_t etaMin=-0.8,
						 Double_t etaMax=0.8,
						 Double_t maxTPCchi2 = -1, 
						 Int_t minNClustersTPC = -1,
						 Bool_t kUsePID = kFALSE,
						 Int_t AODfilterBit = 128,
						 Bool_t bCentralTrigger = kFALSE,
						 TString fileNameBase="AnalysisResults") {

  // Creates a balance function analysis task and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
  TString centralityName("");
  centralityName+=Form("%.0f",centrMin);
  centralityName+="-";
  centralityName+=Form("%.0f",centrMax);
  centralityName+="_";
  centralityName+=Form("%s",centralityEstimator.Data());
  centralityName+="_";
  centralityName+=Form("vZ%.1f",vertexZ);
  centralityName+="_";
  centralityName+=Form("DCAxy%.1f",DCAxy);
  centralityName+="_";
  centralityName+=Form("DCAz%.1f",DCAz);
  centralityName+="_Pt";
  centralityName+=Form("%.1f",ptMin);
  centralityName+="-";
  centralityName+=Form("%.1f",ptMax);
  centralityName+="_Eta";
  centralityName+=Form("%.1f",etaMin);
  centralityName+="-";
  centralityName+=Form("%.1f",etaMax);
  centralityName+="_Chi";
  centralityName+=Form("%.1f",maxTPCchi2);
  centralityName+="_nClus";
  centralityName+=Form("%d",minNClustersTPC);
  centralityName+="_Bit";
  centralityName+=Form("%d",AODfilterBit);
  if(bCentralTrigger)   centralityName+="_withCentralTrigger";



  TString outputFileName(fileNameBase);
  outputFileName.Append(".root");

  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskTriggeredBF", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //===========================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskTriggeredBF", "This task requires an input event handler");
    return NULL;
  }
  TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())) analysisType = "MC";

  // setup the balance function objects
  AliBalanceTriggered *bf  = 0;  // Balance Function object
  AliBalanceTriggered *bfs = 0;  // shuffled Balance function object
  AliBalanceTriggered *bfm = 0;  // mixing Balance function object
  
  if (analysisType=="AOD"){
    
    bf = new AliBalanceTriggered();
    bf->SetAnalysisLevel(analysisType);
    bf->InitHistograms();
    
    if(gRunShuffling){
      bfs = new AliBalanceTriggered();
      bfs->SetAnalysisLevel(analysisType);
      bfs->InitHistograms();
    }
    if(gRunMixing){
      bfm = new AliBalanceTriggered();
      bfm->SetAnalysisLevel(analysisType);
      bfm->InitHistograms();
    }
  }
  else{
    ::Error("AddTaskTriggeredBF", "analysis type NOT supported.");
    return NULL;
  }

  // Create the task, add it to manager and configure it.
  //===========================================================================
  AliAnalysisTaskTriggeredBF *taskTriggeredBF = new AliAnalysisTaskTriggeredBF("TaskTriggeredBF");
  taskTriggeredBF->SetAnalysisObject(bf);
  if(gRunShuffling) taskTriggeredBF->SetShufflingObject(bfs);
  if(gRunMixing) taskTriggeredBF->SetMixingObject(bfm);

  taskTriggeredBF->SetCentralityPercentileRange(centrMin,centrMax);
  if(analysisType == "AOD") {
    // pt and eta cut (pt_min, pt_max, eta_min, eta_max)
    taskTriggeredBF->SetAODtrackCutBit(AODfilterBit);
    taskTriggeredBF->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax);
    
    // set extra DCA cuts (-1 no extra cut)
    taskTriggeredBF->SetExtraDCACutsAOD(DCAxy,DCAz);
    
    // set extra TPC chi2 / nr of clusters cut
    taskTriggeredBF->SetExtraTPCCutsAOD(maxTPCchi2, minNClustersTPC);
    
  }

  // offline trigger selection (AliVEvent.h)
  // taskTriggeredBF->UseOfflineTrigger(); // NOT used (selection is done with the AliAnalysisTaskSE::SelectCollisionCandidates()) 
  // with this only selected events are analyzed (first 2 bins in event QA histogram are the same))
  // documentation in https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PWG1EvSelDocumentation
  if(bCentralTrigger) taskTriggeredBF->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  else                taskTriggeredBF->SelectCollisionCandidates(AliVEvent::kMB);

  // centrality estimator (default = V0M)
  taskTriggeredBF->SetCentralityEstimator(centralityEstimator);
  
  // vertex cut (x,y,z)
  taskTriggeredBF->SetVertexDiamond(.3,.3,vertexZ);
  
  //bf->PrintAnalysisSettings();
  mgr->AddTask(taskTriggeredBF);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCFEbyE.outputBalanceFunctionTriggeredAnalysis";
  AliAnalysisDataContainer *coutQA = mgr->CreateContainer(Form("listQA_%s",centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  AliAnalysisDataContainer *coutTriggeredBF = mgr->CreateContainer(Form("listTriggeredBF_%s",centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  if(gRunShuffling) AliAnalysisDataContainer *coutTriggeredBFS = mgr->CreateContainer(Form("listTriggeredBFShuffled_%s",centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());

  mgr->ConnectInput(taskTriggeredBF, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskTriggeredBF, 1, coutQA);
  mgr->ConnectOutput(taskTriggeredBF, 2, coutTriggeredBF);
  if(gRunShuffling) mgr->ConnectOutput(taskTriggeredBF, 3, coutTriggeredBFS);

  return taskTriggeredBF;
}
