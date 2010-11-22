AliAnalysisTaskMuonFakes* AddTaskMuonFakes(Bool_t useMCLabels = kFALSE)
{
  /// Add AliAnalysisTaskMuonFakes to the train (Philippe Pillot)
  
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTaskMuonFakes","AliAnalysisManager not set!");
    return NULL;
  }
  
  // This task run on ESDs
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD")) {
    Error("AddTaskMuonFakes", "ESD input handler needed!");
    return NULL;
  }
  
  // Create and configure task
  AliAnalysisTaskMuonFakes *task = new AliAnalysisTaskMuonFakes("MUONFakes");
  if (!task) {
    Error("AddTaskMuonFakes", "Muon fakes task cannot be created!");
    return NULL;
  }
  task->UseMCLabels(useMCLabels);
  
  // Add task to analysis manager
  mgr->AddTask(task);
  
  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  // Define output file directory
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  if ( outputfile.IsNull() ) {
    Error("AddTaskMuonFakes", "Common output file is not defined!");
    return NULL;
  }
  outputfile += ":MUON_Fakes";
  
  // Create and connect output containers
  AliAnalysisDataContainer *cout_histo = mgr->CreateContainer("histos", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *cout_track = mgr->CreateContainer("track statistics", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *cout_fakeTrack = mgr->CreateContainer("fake track statistics", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *cout_matchTrack = mgr->CreateContainer("matched track statistics", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *cout_event = mgr->CreateContainer("event statistics", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  mgr->ConnectOutput(task, 1, cout_histo);
  mgr->ConnectOutput(task, 2, cout_track);
  mgr->ConnectOutput(task, 3, cout_fakeTrack);
  mgr->ConnectOutput(task, 4, cout_matchTrack);
  mgr->ConnectOutput(task, 5, cout_event);
  
  return task;
}

