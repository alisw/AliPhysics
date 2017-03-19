// -*- C++ -*-

AliAnalysisTaskDiffCrossSections* AddAnalysisTaskDiffCrossSections(Bool_t isMC,
								   TString mcType,
								   TString triggerSelection,
								   Float_t fmdMultLowCut=0.3)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    mgr = new AliAnalysisManager("My test train");

  if (isMC) {
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }

  AliAnalysisTaskDiffCrossSections::PseudoTrack::Class()->IgnoreTObjectStreamer();
  AliAnalysisTaskDiffCrossSections::PseudoTracks::Class()->IgnoreTObjectStreamer();
  AliAnalysisTaskDiffCrossSections::TreeData::Class()->IgnoreTObjectStreamer();
  AliAnalysisTaskDiffCrossSections::MCInfo::Class()->IgnoreTObjectStreamer();

  AliAnalysisTaskDiffCrossSections* task = new AliAnalysisTaskDiffCrossSections;
  // task->SelectCollisionCandidates(AliVEvent::kMB);

  task->SetIsMC(isMC);
  task->SetMCType(mcType);
  task->SetTriggerSelection(triggerSelection);
  task->SetFMDMultLowCut(fmdMultLowCut);

  AliAnalysisDataContainer* output =
    mgr->CreateContainer(task->GetTreeName(), TTree::Class(), AliAnalysisManager::kOutputContainer,
			 TString(AliAnalysisManager::GetCommonFileName())+":"+task->GetResultsFileName());
  mgr->AddTask(task);

  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);
  return task;
}
