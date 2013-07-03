AliAnalysisTaskEpRatio* AddTaskPHOSEpRatio (Bool_t kMC = kFALSE,
					    const char* name = "PHOSEpRatio",
					    const char* options = "",
					    UInt_t offlineTriggerMask = AliVEvent::kINT7 )
{
  //Add a task AliAnalysisTaskEpRatio to the analysis train.
  //Author: Boris Polishchuk.
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSEpRatio", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSEpRatio", "This task requires an input event handler");
    return NULL;
  }
  
  AliAnalysisTaskEpRatio* task = new AliAnalysisTaskEpRatio(Form("%sTask", name));  
  if(!kMC) task->SelectCollisionCandidates(offlineTriggerMask);
  
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  TString cname(Form("%sCoutput1", name));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
