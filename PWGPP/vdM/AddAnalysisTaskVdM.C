// -*- C++ -*-

AliAnalysisTaskVdM* AddAnalysisTaskVdM(TString branchNames="VertexTracks VertexTracksUnconstrained")
{
  // create manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("VdM analysis manager");

  AliAnalysisTaskVdM* task = new AliAnalysisTaskVdM;

  task->SetBranchNames(branchNames);

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
