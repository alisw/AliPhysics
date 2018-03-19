AliAnalysisTaskPHOSTriggerQA* AddTaskPHOSTriggerQA(TString fname="PHOSTriggerQA.root", TString contname="")
{
  //Add PHOS trigger QA task to the PWGPP QA train.

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSTriggerQA", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSTriggerQA", "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisTaskPHOSTriggerQA* task = new AliAnalysisTaskPHOSTriggerQA("PHOSTriggerQA");
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  // container output into particular file
  if (fname.Length() != 0 && contname.Length() != 0)
    mgr->ConnectOutput(task, 1, mgr->CreateContainer(contname.Data(),TList::Class(), AliAnalysisManager::kOutputContainer, fname.Data()));
  
  // container output into common file
  if (fname.Length() == 0) {
    if (contname.Length() == 0) contname = "PHOSTriggerQAResults";
    mgr->ConnectOutput(task, 1, mgr->CreateContainer(contname.Data(),TList::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName()));      
  }
  
  return task;
}
