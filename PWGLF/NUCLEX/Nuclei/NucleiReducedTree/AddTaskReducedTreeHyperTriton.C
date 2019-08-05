AliAnalysisTask *AddTaskReducedTreeHyperTriton ()  {

  // options central SemiCentral and kINT7
  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("No analysis manager found.");
    return 0;
  }

  AliAnalysisTaskReducedTreeHypertriton *task = new AliAnalysisTaskReducedTreeHypertriton("TaskHypTri");
  
  // Min bais
  task -> SelectCollisionCandidates (AliVEvent::kINT7);
  task -> SetCentrality (0.0,90.0);

//   // central
//   task -> SelectCollisionCandidates (AliVEvent::kCentral);
//   task -> SetCentrality (0.0,10.0);
//     
//   // SemiCentral Collision
//   task -> SelectCollisionCandidates (AliVEvent::kSemiCentral);
//   task -> SetCentrality (30.0,50.0);
  
  
  mgr -> AddTask(task);
  TString Filename =mgr->GetCommonFileName();
  AliAnalysisDataContainer *cQA = mgr->CreateContainer("QAHistograms", TList::Class(), AliAnalysisManager::kOutputContainer, Filename.Data());
  AliAnalysisDataContainer *cOutput = mgr->CreateContainer("ReducedTree_Hypertriton", TList::Class(), AliAnalysisManager::kOutputContainer, Filename.Data());

  mgr -> ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  mgr -> ConnectOutput (task, 1, cOutput);
  mgr -> ConnectOutput (task, 2, cQA);
    
  return task;
    
}