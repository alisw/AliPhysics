AliAnalysisTaskMuonFakes* AddTaskMuonFakes(Bool_t useMCLabels = kFALSE, Bool_t combineMCId = kFALSE,
					   Bool_t matchTrig = kFALSE, Bool_t applyAccCut = kFALSE,
					   TString extension = "")
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
  TString name = Form("MUONFakes%s",extension.Data());
  AliAnalysisTaskMuonFakes *task = new AliAnalysisTaskMuonFakes(name.Data());
  if (!task) {
    Error("AddTaskMuonFakes", "Muon fakes task cannot be created!");
    return NULL;
  }
  task->UseMCLabels(useMCLabels);
  task->CombineMCId(combineMCId);
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
  TString suffix = (!extension.IsNull()) ? Form("_%s",extension.Data()) : "";
  outputfile += suffix;
  
  // Create and connect output containers
  AliAnalysisDataContainer *cout_histo = mgr->CreateContainer(Form("histos%s",suffix.Data()), TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *cout_track = mgr->CreateContainer(Form("trackCounters%s",suffix.Data()), AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *cout_fakeTrack = mgr->CreateContainer(Form("fakeTrackCounters%s",suffix.Data()), AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *cout_matchTrack = mgr->CreateContainer(Form("matchedTrackCounters%s",suffix.Data()), AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *cout_event = mgr->CreateContainer(Form("eventCounters%s",suffix.Data()), AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *cout_histo2 = mgr->CreateContainer(Form("histos2%s",suffix.Data()), TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *cout_pair = mgr->CreateContainer(Form("pairCounters%s",suffix.Data()), AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
  mgr->ConnectOutput(task, 1, cout_histo);
  mgr->ConnectOutput(task, 2, cout_track);
  mgr->ConnectOutput(task, 3, cout_fakeTrack);
  mgr->ConnectOutput(task, 4, cout_matchTrack);
  mgr->ConnectOutput(task, 5, cout_event);
  mgr->ConnectOutput(task, 6, cout_histo2);
  mgr->ConnectOutput(task, 7, cout_pair);
  
  return task;
}

