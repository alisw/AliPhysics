AliAnalysisTaskPHOSPbPbQARun2* AddTaskPHOSPbPbQARun2( Bool_t isMC = kFALSE)
{
  //Add PHOS PbPb QA task to the PWGPP QA train.
  //See PHOSPbPb.C how to run it locally or standalone.

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSPbPbQA", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSPbPbQA", "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisTaskPHOSPbPbQARun2* task = new AliAnalysisTaskPHOSPbPbQARun2("PbPbQA");
  if (!task) return NULL;

  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("PbPbQA", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:PbPbQA",AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectOutput(task, 1, coutput1);
  return task;
}
