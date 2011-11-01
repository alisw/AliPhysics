AliAnalysisTaskPHOSPbPbQA* AddTaskPHOSPbPbQA(char* fname="PHOSPbPbQA.root",
					     char* contname=NULL)
{
  //Add PHOS PbPb QA task to the PWG1 QA train.
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

  AliAnalysisTaskPHOSPbPbQA* task = new AliAnalysisTaskPHOSPbPbQA("PbPbQA");
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  // container output into particular file
  if (fname && contname)
    mgr->ConnectOutput(task, 1, mgr->CreateContainer(contname,TList::Class(), AliAnalysisManager::kOutputContainer, fname));
  
  // container output into common file
  if (!fname) {
    if (!contname) contname = "PHOSPbPbQAResults";
    mgr->ConnectOutput(task, 1, mgr->CreateContainer(contname,TList::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName()));		       
  }
  
  return task;
}
