AliAnalysisTaskGenTuner* AddTaskGenTuner()
{
  /// Add AliAnalysisTaskGenTuner to the train (Philippe Pillot)
  
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTaskGenTuner","AliAnalysisManager not set!");
    return NULL;
  }
  
  // This task run on ESDs or AODs
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("AOD")) {
    Error("AddTaskGenTuner", "ESD input handler needed!");
    return NULL;
  }
  
  // Create and configure task
  AliAnalysisTaskGenTuner *task = new AliAnalysisTaskGenTuner("GenTuner");
  if (!task) {
    Error("AddTaskGenTuner", "Muon physics task cannot be created!");
    return NULL;
  }
  
  // Add task to analysis manager
  mgr->AddTask(task);
  
  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  // Define output file directory
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  if ( outputfile.IsNull() ) {
    Error("AddTaskGenTuner", "Common output file is not defined!");
    return NULL;
  }
  outputfile += ":MUON_GenTuner";
  
  // Create and connect output containers
  AliAnalysisDataContainer *histo = mgr->CreateContainer("Histograms", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *eventStat = mgr->CreateContainer("eventCounters", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  mgr->ConnectOutput(task, 1, histo);
  mgr->ConnectOutput(task, 2, eventStat);
  
  return task;
}

