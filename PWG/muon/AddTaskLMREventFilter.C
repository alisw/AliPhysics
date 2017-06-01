
AliAnalysisTask *AddTaskLMREventFilter()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      ::Error("AddTaskLMRMuMu", "No analysis manager to connect to.");
      return NULL;
    }
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskLMRMuMu", "This task requires an input event handler");
      return NULL;
    }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  

  AliMuonTrackCuts* MuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts", "StandardMuonTrackCuts"); 
  MuonTrackCuts->Print("all");

  AliAnalysisTaskLMREventFilter* task = new AliAnalysisTaskLMREventFilter("AliAnalysisTaskLMREventFilter", MuonTrackCuts);

  mgr->AddTask(task);

  AliAnalysisDataContainer *output1 = mgr->CreateContainer("list"    , TList::Class(), AliAnalysisManager::kOutputContainer, "AliLMREvents.root");
  AliAnalysisDataContainer *output2 = mgr->CreateContainer("MuonTree", TTree::Class(), AliAnalysisManager::kOutputContainer, "AliLMREvents.root");
  
  // finaly connect input and output
  mgr->ConnectInput(task, 0,  mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output1);
  mgr->ConnectOutput(task, 2, output2);

  return task;
}
