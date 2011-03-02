AliAnalysisTaskITSAlignQA *AddTaskITSAlign(Int_t nrun, Int_t year){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskITSAlign", "No analysis manager to connect to.");
    return NULL;
  }   
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskITSAlign", "This task requires an input event handler");
    return NULL;
  }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); 
  if(type.Contains("AOD")){
    ::Error("AddTaskITSAlign", "This task requires to run on ESD");
    return NULL;
  }
  
  // Create and configure the task
  AliAnalysisTaskITSAlignQA *taskali = new AliAnalysisTaskITSAlignQA();
  //  taskali->SelectCollisionCandidates();
  taskali->SetOCDBInfo(nrun,Form("alien://folder=/alice/data/%d/OCDB",year)) ; 
  mgr->AddTask(taskali);
  
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":ITSAlignQA";
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clistITSAlignQA",
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outputFileName );
  
  mgr->ConnectInput(taskali, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskali, 1, coutput1);
  return taskali;
}
