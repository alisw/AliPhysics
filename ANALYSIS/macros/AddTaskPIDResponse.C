AliAnalysisTask *AddTaskPIDResponse()
{
// Macro to connect a centrality selection task to an existing analysis manager.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPIDResponse", "No analysis manager to connect to.");
    return 0x0;
  }

  Bool_t isMC=kFALSE;
  if (mgr->GetInputEventHandler()->IsA() == AliESDInputHandler::Class()) {
    isMC=mgr->GetMCtruthEventHandler()!=0x0;
  }
  
  AliAnalysisTaskPIDResponse *pidTask = new AliAnalysisTaskPIDResponse("PIDResponseTask");
//   pidTask->SelectCollisionCandidates(AliVEvent::kMB);
  pidTask->SetIsMC(isMC);
  mgr->AddTask(pidTask);
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("PIDResponseQA",
    TList::Class(), AliAnalysisManager::kOutputContainer,
    "PIDResponseQA.root");
  
  mgr->ConnectInput(pidTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(pidTask,1,coutput1);
  
  return pidTask;
}
