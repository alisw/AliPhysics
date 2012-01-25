AliAnalysisTaskMuonFakes* AddTaskMuonFakes(Bool_t useMCLabels = kFALSE, Bool_t matchTrig = kFALSE, Bool_t applyAccCut = kFALSE)
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
  task->MatchTrigger(matchTrig);
  task->ApplyAccCut(applyAccCut);
  
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
  AliAnalysisDataContainer *cout_track = mgr->CreateContainer("trackCounters", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *cout_fakeTrack = mgr->CreateContainer("fakeTrackCounters", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *cout_matchTrack = mgr->CreateContainer("matchedTrackCounters", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *cout_event = mgr->CreateContainer("eventCounters", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *cout_histo2 = mgr->CreateContainer("histos2", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *cout_pair = mgr->CreateContainer("pairCounters", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  mgr->ConnectOutput(task, 1, cout_histo);
  mgr->ConnectOutput(task, 2, cout_track);
  mgr->ConnectOutput(task, 3, cout_fakeTrack);
  mgr->ConnectOutput(task, 4, cout_matchTrack);
  mgr->ConnectOutput(task, 5, cout_event);
  mgr->ConnectOutput(task, 6, cout_histo2);
  mgr->ConnectOutput(task, 7, cout_pair);
  
  return task;
}

