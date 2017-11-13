AliTaskMuonTrackSmearingQA* AddTaskMuonTrackSmearingQA(Int_t chosenFunction = 0)
{
  /// Add AliTaskMuonTrackSmearingQA to the train (Philippe Pillot)
  
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    Error("AddTaskMuonTrackSmearing","AliAnalysisManager not set!");
    return NULL;
  }
  
  // This task run on ESDs or AODs
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    Error("AddTaskMuonTrackSmearing", "ESD or AOD input handler needed!");
    return NULL;
  }
  
  // Create and configure task
  AliTaskMuonTrackSmearingQA *task = new AliTaskMuonTrackSmearingQA("MuonTrackSmearingQA",chosenFunction);
  if (!task) {
    Error("AddTaskMuonTrackSmearing", "Track smearing QA task cannot be created!");
    return NULL;
  }

  // Add task to analysis manager
  mgr->AddTask(task);
  
  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  // Define output file directory
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  if ( outputfile.IsNull() ) {
    Error("AddTaskMuonTrackSmearing", "Common output file is not defined!");
    return NULL;
  }
  outputfile += ":MUON_TrackSmearingQA";
  
  // Create and connect output containers
  AliAnalysisDataContainer *genHistos = mgr->CreateContainer("GenHistos", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *recHistos = mgr->CreateContainer("RecHistos", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *resHistos = mgr->CreateContainer("ResHistos", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  mgr->ConnectOutput(task, 1, genHistos);
  mgr->ConnectOutput(task, 2, recHistos);
  mgr->ConnectOutput(task, 3, resHistos);

   return task;
}
