AliAnalysisTaskMuonQA *AddTaskMuonQA(Bool_t selectPhysics = kTRUE, Bool_t selectMatched = kTRUE, 
				     Bool_t applyAccCut = kTRUE, Bool_t selectTrigger = kFALSE,
				     UInt_t triggerMask = AliVEvent::kMUS7,
						 Short_t selectCharge = 0)
{
  /// Add AliAnalysisTaskMuonQA to the train (Philippe Pillot)
  
  
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTaskMuonQA","AliAnalysisManager not set!");
    return NULL;
  }
  
  // This task run on ESDs
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD")) {
    Error("AddTaskMuonQA", "ESD input handler needed!");
    return NULL;
  }
  
  // Create and configure task
  AliAnalysisTaskMuonQA *task = new AliAnalysisTaskMuonQA("MuonQA");
  if (!task) {
    Error("AddTaskMuonQA", "Muon QA task cannot be created!");
    return NULL;
  }
  task->SelectPhysics(selectPhysics);
  task->SelectTrigger(selectTrigger, triggerMask);
  task->SelectMatched(selectMatched);
  task->ApplyAccCut(applyAccCut);
  task->SelectCharge(selectCharge);
  
  // Add task to analysis manager
  mgr->AddTask(task);
  
  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  // Define output file directory
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  if ( outputfile.IsNull() ) {
    Error("AddTaskMuonQA", "Common output file is not defined!");
    return NULL;
  }
  outputfile += ":MUON_QA";
  
  // Create and connect output containers
  AliAnalysisDataContainer *cout_histo1 = mgr->CreateContainer("general1", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *cout_histo2 = mgr->CreateContainer("expert", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *cout_trackStat = mgr->CreateContainer("trackCounters", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *cout_eventStat = mgr->CreateContainer("eventCounters", AliCounterCollection::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *cout_normalized = mgr->CreateContainer("general2", TObjArray::Class(), AliAnalysisManager::kParamContainer, outputfile);
  mgr->ConnectOutput(task, 1, cout_histo1);
  mgr->ConnectOutput(task, 2, cout_histo2);
  mgr->ConnectOutput(task, 3, cout_trackStat);
  mgr->ConnectOutput(task, 4, cout_eventStat);
  mgr->ConnectOutput(task, 5, cout_normalized);
  
  return task;
}   
