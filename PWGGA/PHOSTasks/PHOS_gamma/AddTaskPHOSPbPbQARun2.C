AliAnalysisTaskPHOSPbPbQARun2* AddTaskPHOSPbPbQARun2(const char* fname="AnalysisResults_PHOSPbPbQA.root",
					             const char* contname=NULL)
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
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  if(!contname) contname = "PHOSPbPbQAResults";
  if(!fname)    fname =  mgr->GetCommonFileName();

  mgr->ConnectOutput(task, 1, mgr->CreateContainer(contname,TList::Class(), AliAnalysisManager::kOutputContainer, fname));
  return task;
  
}
