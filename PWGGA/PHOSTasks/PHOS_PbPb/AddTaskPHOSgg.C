AliAnalysisTaskgg* AddTaskPHOSgg(const char* name = "PHOSgg",
                                 UInt_t offlineTriggerMask = AliVEvent::kCentral,
				 Bool_t isPbPb=kTRUE)
{
  //Add a task AliAnalysisTaskGammaFlow to the analysis train
  //Author: Dmitri Peressounko
  /* $Id$ */

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSgg", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSgg", "This task requires an input event handler");
    return NULL;
  }
 
  AliAnalysisTaskgg* task = new AliAnalysisTaskgg(Form("%sTask", name));
  task->SelectCollisionCandidates(offlineTriggerMask);
  task->SelectPbPb(isPbPb) ;
  
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer(); 
  TString trigger;
  if(offlineTriggerMask==AliVEvent::kCentral) trigger = "Central" ;
  else if(offlineTriggerMask==AliVEvent::kSemiCentral) trigger = "SemiCentral" ;
  else trigger = "Other" ;
  TString cname(Form("%s%s", name,trigger.Data()));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectInput(task , 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
