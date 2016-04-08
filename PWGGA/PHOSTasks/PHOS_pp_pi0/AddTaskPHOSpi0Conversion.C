void AddTaskPHOSpi0Conversion(const char* name = "PHOSconv")
{
  //Add a task AliAnalysisTaskPi0Conversion to the analysis train
  //Author: Yulia Demchenko
  /* $Id$ */

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSGammaFlow", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSGammaFlow", "This task requires an input event handler");
    return NULL;
  }
  
  AliAnalysisTaskPi0Conversion * task = new AliAnalysisTaskPi0Conversion(Form("%s", name));
  
  task->SelectCollisionCandidates(AliVEvent::kMB);
  
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer(); 
  TString cname(Form("%s", name));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectInput(task , 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
