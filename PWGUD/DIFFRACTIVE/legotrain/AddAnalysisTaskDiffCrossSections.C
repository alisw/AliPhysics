// -*- C++ -*-

AliAnalysisTaskDiffCrossSections* AddAnalysisTaskDiffCrossSections(Bool_t isMC,
								   TString mcType,
								   TString triggerSelection)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    mgr = new AliAnalysisManager("My test train");
  
  if (isMC) {
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }

  AliAnalysisTaskDiffCrossSections* task = new AliAnalysisTaskDiffCrossSections;
  // task->SelectCollisionCandidates(AliVEvent::kMB);

  task->SetIsMC(isMC);
  task->SetMCType(mcType);
  task->SetTriggerSelection(triggerSelection);

  AliAnalysisDataContainer* output =
    mgr->CreateContainer(task->GetTreeName(), TTree::Class(), AliAnalysisManager::kOutputContainer,
			 TString(AliAnalysisManager::GetCommonFileName())+":"+task->GetResultsFileName());  
  mgr->AddTask(task);
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);
  return task;
}

