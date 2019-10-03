AliAnalysisTaskNorm *AddTaskNorm ( Bool_t isESD = kFALSE, Bool_t isMC = kFALSE, TString sBeamConf = "p-Pb", TString sYear = "2013" )
{
  /// Add AliAnalysisTaskNorm to the train 
  Info("AddTaskNorm", Form("beamConf %s isESD %d, isMC %d", sBeamConf.Data(), isESD, isMC));
  
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
    Error("AddTaskNorm","AliAnalysisManager not set!");
    return NULL;
  }
  
  // This task run on AOD and ESD
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("AOD") && !isESD) {
    Error("AddTaskNorm", "AOD input handler needed!");
    return NULL;
  }
  if ( !type.Contains("ESD") && isESD ) {
    Error("AddTaskNorm","ESD input handler needed!");
    return NULL;
  }

  AliAnalysisTaskNorm *task = new AliAnalysisTaskNorm("TaskNorm");
  if (!task) {
    Error("AddTaskNorm", "Analysis task Norm cannot be created!");
    return NULL;
  }
  task->SetIsMC(isMC);
  task->SetIsESD(isESD);
  task->SetBeamConf(sBeamConf);
  ((AliMuonEventCuts*) task->GetMuonEventCuts())->SetFilterMask ( AliMuonEventCuts::kSelectedTrig );
  // 2013 p-Pb and Pb-p
  if ( sYear.Contains("2013") )   ((AliMuonEventCuts*) task->GetMuonEventCuts())->SetTrigClassPatterns("CINT7-B-NOPF-ALLNOTRD,CINT7-B-NOPF-ALLNOTRD&0MSL,CINT7-B-NOPF-ALLNOTRD&0MUL,CMSL7-B-NOPF-MUON,CMSL7-B-NOPF-MUON&0MUL,CMUL7-B-NOPF-MUON,CMUL7-B-NOPF-ALLNOTRD","0MSL:12,0MUL:14");
  //2015 Pb-Pb
  if ( sYear.Contains("2015") ) ((AliMuonEventCuts*) task->GetMuonEventCuts())->SetTrigClassPatterns("CINT7-B-NOPF-MUFAST,CINT7-B-NOPF-MUFAST&0MSL,CINT7-B-NOPF-MUFAST&0MUL,CMSL7-B-NOPF-MUFAST,CMSL7-B-NOPF-MUFAST&0MUL,CMUL7-B-NOPF-MUFAST","0MSL:12,0MUL:14");
  // Add task to analysis manager
  mgr->AddTask(task);
  
  // Connect input container
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  // Define output file directory
  TString outputFile = AliAnalysisManager::GetCommonFileName();
  if ( outputFile.IsNull() ) {
    Error("AddTaskMuon", "Common output file is not defined!");
    return NULL;
  }

  if ( isESD )
    outputFile = "Norm.ESD.Rec.root";
  else 
    outputFile = "Norm.AOD.Rec.root";

  // Create and connect output containers
  AliAnalysisDataContainer *coutputeventcounter = mgr->CreateContainer("eventCounters", AliCounterCollection::Class(),    AliAnalysisManager::kOutputContainer, outputFile);
  AliAnalysisDataContainer *coutputruncounter = mgr->CreateContainer("runCounters", AliCounterCollection::Class(),    AliAnalysisManager::kOutputContainer, outputFile);
  AliAnalysisDataContainer *coutputvertexlist = mgr->CreateContainer("vertex", TObjArray::Class(),    AliAnalysisManager::kOutputContainer, outputFile);
  AliAnalysisDataContainer *coutputv0alist = mgr->CreateContainer("V0A", TObjArray::Class(),    AliAnalysisManager::kOutputContainer, outputFile);
  AliAnalysisDataContainer *coutputznlist = mgr->CreateContainer("ZN", TObjArray::Class(),    AliAnalysisManager::kOutputContainer, outputFile);
  
  mgr->ConnectOutput(task, 1, coutputeventcounter);
  mgr->ConnectOutput(task, 2, coutputruncounter);
  mgr->ConnectOutput(task, 3, coutputvertexlist);
  mgr->ConnectOutput(task, 4, coutputv0alist);
  mgr->ConnectOutput(task, 5, coutputznlist);

  return task;
}   
