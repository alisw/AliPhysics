AliAnalysisTaskSimpleTreeMaker *AddTaskSimpleTreeMaker(TString taskName = "MLtree", 
                                             Double_t etaMin = -0.8,
                                             Double_t etaMax = 0.8,
                                             Double_t ptMin = 0.2,
                                             Double_t ptMax = 10.0,
					     ) {				

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskSimpleTreeMaker",  "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //===========================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskSimpleTreeMaker",  "This task requires an input event handler");
    return NULL;
  }

  TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  if (analysisType!="ESD"){
    ::Error("AddTaskSimpleTreeMaker",  "analysis type NOT AOD --> makes no sense!");
    return NULL;
  }

  AliAnalysisTaskSimpleTreeMaker *taskESD = new AliAnalysisTaskSimpleTreeMaker(taskName);
  // ==========================================================================
  // user customization part

  taskESD->SelectCollisionCandidates(AliVEvent::kINT7);
  //taskESD->SetMC(kFALSE);
  //taskESD->setSDDstatus(kFALSE);
  //taskESD->createV0tree(kFALSE);


  
  // ==========================================================================
  mgr->AddTask(taskESD);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
 
  TString outputFileName = AliAnalysisManager::GetCommonFileName();  

  AliAnalysisDataContainer *coutTree = mgr->CreateContainer("Tree", TTree::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  AliAnalysisDataContainer *coutHisto1 = mgr->CreateContainer("Histo", TH1F::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  AliAnalysisDataContainer *coutHisto2 = mgr->CreateContainer("Arm. Plot", TH2F::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());

  mgr->ConnectInput(taskESD, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskESD, 1, coutTree);
  mgr->ConnectOutput(taskESD, 2, coutHisto1);
  mgr->ConnectOutput(taskESD, 3, coutHisto2);

 return taskESD;
}
