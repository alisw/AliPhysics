class AliAnalysisDataContainer;
AliAnalysisTaskEtaPhigg* AddTaskPHOSEtaPhi(const char* name = "PHOSEtaPhi",
                                           UInt_t offlineTriggerMask = AliVEvent::kCentral)
{
  AliAnalysisTaskEtaPhigg * task = new AliAnalysisTaskEtaPhigg(Form("%sTask", name));
  task->SelectCollisionCandidates(offlineTriggerMask);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer(); 
  TString cname(Form("%s", name));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectInput(task , 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
