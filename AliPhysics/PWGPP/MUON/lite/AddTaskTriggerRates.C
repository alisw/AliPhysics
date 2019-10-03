AliAnalysisTaskTriggerRates* AddTaskTriggerRates(TString extension = "")
{
  /// Add AliAnalysisTaskTriggerRates to the train (Philippe Pillot)
  
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    Error("AddTaskTriggerRates","AliAnalysisManager not set!");
    return NULL;
  }
  
  // This task run on ESDs
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD")) {
    Error("AddTaskTriggerRates", "ESD input handler needed!");
    return NULL;
  }
  
  // Create and configure task
  TString suffix = (!extension.IsNull()) ? Form("_%s",extension.Data()) : "";
  AliAnalysisTaskTriggerRates *task = new AliAnalysisTaskTriggerRates(Form("TriggerRates%s",suffix.Data()));
  if (!task) {
    Error("AddTaskTriggerRates", "Trigger rates task cannot be created!");
    return NULL;
  }
  
  // Add task to analysis manager
  mgr->AddTask(task);
  
  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  // Define output file directory
  TString outputfile = "Output.root"; // mgr->GetCommonFileName();
  if ( outputfile.IsNull() ) {
    Error("AddTaskTriggerRates", "Common output file is not defined!");
    return NULL;
  }
  
  // Create and connect output containers
  AliAnalysisDataContainer *trigStat = mgr->CreateContainer(Form("triggerCounters%s",suffix.Data()), AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  mgr->ConnectOutput(task, 1, trigStat);
  
  return task;
}

