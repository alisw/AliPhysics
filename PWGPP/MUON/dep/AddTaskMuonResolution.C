AliAnalysisTaskMuonResolution *AddTaskMuonResolution(Bool_t selectPhysics = kFALSE, Bool_t selectTrigger = kFALSE,
						     Bool_t matchTrig = kTRUE, Bool_t applyAccCut = kTRUE,
						     Double_t minMomentum = 0., Bool_t correctForSystematics = kTRUE,
						     Int_t extrapMode = 1)
{
  /// Add AliAnalysisTaskMuonResolution to the train (Philippe Pillot)
  
  
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTaskMuonResolution","AliAnalysisManager not set!");
    return NULL;
  }
  
  // This task run on ESDs
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD")) {
    Error("AddTaskMuonResolution", "ESD input handler needed!");
    return NULL;
  }
  
  // Create and configure task
  AliAnalysisTaskMuonResolution *task = new AliAnalysisTaskMuonResolution("MuonResolution");
  if (!task) {
    Error("AddTaskMuonResolution", "Muon resolution task cannot be created!");
    return NULL;
  }
  task->SelectPhysics(selectPhysics);
  task->SelectTrigger(selectTrigger);
//  task->SelectTrigger(selectTrigger, AliVEvent::kMB);
  task->MatchTrigger(matchTrig);
  task->ApplyAccCut(applyAccCut);
  task->SetMinMomentum(minMomentum);
  task->CorrectForSystematics(correctForSystematics);
  task->SetExtrapMode(extrapMode);
  
  // Add task to analysis manager
  mgr->AddTask(task);
  
  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  // Define output file directory
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  if ( outputfile.IsNull() ) {
    Error("AddTaskMuonResolution", "Common output file is not defined!");
    return NULL;
  }
  outputfile += ":MUON_Resolution";
  
  // Create and connect output containers
  AliAnalysisDataContainer *cout_histo1 = mgr->CreateContainer("Residuals", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *cout_histo2 = mgr->CreateContainer("ResidualsVsP", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *cout_histo3 = mgr->CreateContainer("ResidualsVsCent", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  AliAnalysisDataContainer *cout_histo4 = mgr->CreateContainer("LocalChi2", TObjArray::Class(), AliAnalysisManager::kParamContainer, outputfile);
  AliAnalysisDataContainer *cout_histo5 = mgr->CreateContainer("ChamberRes", TObjArray::Class(), AliAnalysisManager::kParamContainer, outputfile);
  AliAnalysisDataContainer *cout_histo6 = mgr->CreateContainer("TrackRes", TObjArray::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  mgr->ConnectOutput(task, 1, cout_histo1);
  mgr->ConnectOutput(task, 2, cout_histo2);
  mgr->ConnectOutput(task, 3, cout_histo3);
  mgr->ConnectOutput(task, 4, cout_histo4);
  mgr->ConnectOutput(task, 5, cout_histo5);
  mgr->ConnectOutput(task, 6, cout_histo6);
  
  return task;
}   
