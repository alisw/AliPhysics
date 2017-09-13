AliAnalysisTaskCheckVertexAOD *AddTaskCheckVertexAOD(TString suffix="",
						     Double_t maxMult=500.,
						     Int_t mask=AliVEvent::kAnyINT,
						     Bool_t readMC=kFALSE)
{

  // Creates, configures and attaches to the train the task for AOD vertex checks
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCheckVertexAOD", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskCheckVertexAOD", "This task requires an input event handler");
    return NULL;
  }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("ESD")){
    ::Error("AddTaskCheckVertexAOD", "This task requires to run on AOD");
    return NULL;
  }
  
  //Bool_t isMC=kFALSE;
  //if (mgr->GetMCtruthEventHandler()) isMC=kTRUE;
  
  // Add MC handler (for kinematics)
  if(readMC){
    AliMCEventHandler* handler = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
    if (!handler) {
      ::Error("AddTaskCheckVertexAOD","Macro called with readMC=true but MC handler not present");
      return 0;
    }
  }
  // Create and configure the task
  AliAnalysisTaskCheckVertexAOD *task = new AliAnalysisTaskCheckVertexAOD();
  task->SetReadMC(readMC);
  task->SetUpperMultiplicity(maxMult);
  if(mask>=0){
    task->SetUsePhysicsSelection(kTRUE);
    task->SetTriggerMask(mask);
  }else{
    task->SetUsePhysicsSelection(kFALSE);
  }
  mgr->AddTask(task);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":CheckVertexAOD";

  TString listname="clistCheckVertexAOD";
  listname+=suffix.Data();
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(listname,
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFileName.Data() );

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  return task;
}   

