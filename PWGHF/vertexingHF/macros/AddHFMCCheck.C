AliAnalysisTaskCheckHFMCProd *AddHFMCCheck(Bool_t isPbPb=kTRUE, Bool_t readMC=kTRUE){

  // Creates, configures and attaches to the train the task for QA of ITS standalone tracks
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AliAnalysisTaskCheckHFMCProd", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AliAnalysisTaskCheckHFMCProd", "This task requires an input event handler");
    return NULL;
  }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD")){
    ::Error("AliAnalysisTaskCheckHFMCProd", "This task requires to run on ESD");
    return NULL;
  }
  
  //Bool_t isMC=kFALSE;
  //if (mgr->GetMCtruthEventHandler()) isMC=kTRUE;
  
  // Add MC handler (for kinematics)
  if(readMC){
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }
  // Create and configure the task
  AliAnalysisTaskCheckHFMCProd *task = new AliAnalysisTaskCheckHFMCProd();
  if(isPbPb) task->SetPbPb();
  task->SetReadMC(readMC);
  mgr->AddTask(task);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":HFMCCheck";
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clistHFMCCheck",
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFileName );
  
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  return task;
}   
