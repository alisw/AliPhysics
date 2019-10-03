// -*- C++ -*-

AliAnalysisTaskDG* AddAnalysisTaskDG(Bool_t isMC,
				     TString branchNames,
				     TString trackCutType,
				     TString triggerSelection,
				     TString triggerSelectionSPD="",
				     Int_t   maxTracksSave=4)
{
  // create manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("My test train");

  if (isMC) {
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }

  AliAnalysisTaskDG* task = new AliAnalysisTaskDG;
  // task->SelectCollisionCandidates(AliVEvent::kMB);

  task->SetIsMC(isMC);
  task->SetBranchNames(branchNames);
  task->SetTrackCutType(trackCutType);
  task->SetTriggerSelection(triggerSelection);
  task->SetTriggerSelectionSPD(triggerSelectionSPD);
  task->SetMaxTracksSave(maxTracksSave);

  Printf("created task");

  // OUTPUT --------------------------------------------------------------------
  AliAnalysisDataContainer* output1 =
    mgr->CreateContainer(task->GetListName(), TList::Class(), AliAnalysisManager::kOutputContainer,
			 TString(AliAnalysisManager::GetCommonFileName())+":"+task->GetResultsFileName());
  AliAnalysisDataContainer* output2 =
    mgr->CreateContainer(task->GetTreeName(), TTree::Class(), AliAnalysisManager::kOutputContainer,
			 TString(AliAnalysisManager::GetCommonFileName())+":"+task->GetResultsFileName());
  mgr->AddTask(task);

  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output1);
  mgr->ConnectOutput(task, 2, output2);
  Printf("set up task connections");
  Printf("--------------------------------------------------------------------------------");
  return task;
}
