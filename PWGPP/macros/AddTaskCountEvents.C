AliAnalysisTaskCountEvents *AddTaskCountEvents(TString suffix="")
{

  // Creates, configures and attaches to the train the task to count events
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCountEvents", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskCountEvents", "This task requires an input event handler");
    return NULL;
  }   
  
  // Create and configure the task
  AliAnalysisTaskCountEvents *task = new AliAnalysisTaskCountEvents();
  mgr->AddTask(task);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":CountEvents";

  TString listname="clistCountEvents";
  listname+=suffix.Data();
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(listname,
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFileName.Data() );

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  return task;
}   

