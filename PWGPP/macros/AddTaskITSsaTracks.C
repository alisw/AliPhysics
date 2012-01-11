AliAnalysisTaskITSsaTracks *AddTaskITSsaTracks(Bool_t readMC=kFALSE,Bool_t UseMCtruthForPID=kFALSE){
  // Creates, configures and attaches to the train the task for pi, K , p spectra
  // with ITS standalone tracks
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskITSsaTracks", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskITSsaTracks", "This task requires an input event handler");
    return NULL;
  }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD")){
    ::Error("AddTaskITSsaTracks", "This task requires to run on ESD");
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
  AliAnalysisTaskITSsaTracks *taskits = new AliAnalysisTaskITSsaTracks();
  taskits->SelectCollisionCandidates();
  taskits->SetMinITSPoints(4);
  taskits->SetReadMC(readMC);
  taskits->SetUseMCtruthForPID(UseMCtruthForPID);
  mgr->AddTask(taskits);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":ITSsaTracks";
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clistITSsaTracks",
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFileName );
  
  mgr->ConnectInput(taskits, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskits, 1, coutput1);
  return taskits;
}   

