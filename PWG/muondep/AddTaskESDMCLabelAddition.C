AliAnalysisTaskESDMCLabelAddition *AddTaskESDMCLabelAddition()
{
  /// Add AliAnalysisTaskESDMCLabelAddition to the train (Philippe Pillot)
  
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTaskESDMCLabelAddition","AliAnalysisManager not set!");
    return NULL;
  }
  
  // This task runs on ESDs
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD")) {
    Error("AddTaskESDMCLabelAddition", "ESD input handler needed!");
    return NULL;
  }
  
  // This task needs MC input
  AliMCEventHandler *mcH = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
  if (!mcH) {
    Error("AddTaskESDMCLabelAddition", "No MC handler connected!");
    return NULL;
  }   
  
  // Create and configure task
  AliAnalysisTaskESDMCLabelAddition *task = new AliAnalysisTaskESDMCLabelAddition("ESD MC Labels addition");
  if (!task) {
    Error("AddTaskESDMCLabelAddition", "MClabel addition task cannot be created!");
    return NULL;
  }
  
  // Add task to analysis manager
  mgr->AddTask(task);
  
  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  return task;
  
}

