AliAnalysisTask *AddTaskReducedTreeNuclei_mh(bool doTri=false, bool doHypTri=false)  {
   //Get the current analysis manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("No analysis manager found.");
      return 0;
   }
   
   AliAnalysisTaskReducedTreeNuclei *task = new AliAnalysisTaskReducedTreeNuclei("TaskNuclei");
//    task->SelectCollisionCandidates(AliVEvent::kINT7/* + AliVEvent::kHighMultV0*/);
// //    task->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kHighMultV0);
   
      task->SelectCollisionCandidates(AliVEvent::kAny);
   
   
   task->FillTritonTree(doTri);
   task->FillHypTritonTree(doHypTri);
   
   mgr->AddTask(task);
   TString Filename = mgr->GetCommonFileName();
   AliAnalysisDataContainer *cQA     = mgr->CreateContainer("QAHistograms", TList::Class(),AliAnalysisManager::kOutputContainer,Filename.Data());
   AliAnalysisDataContainer *cOutput = mgr->CreateContainer("Results"     , TList::Class(),AliAnalysisManager::kOutputContainer,Filename.Data());
   
   Filename	+= ":Trees";
   AliAnalysisDataContainer *cHeliumTree      = mgr->CreateContainer("reducedTree_Helium",      TTree::Class(),AliAnalysisManager::kOutputContainer,Filename.Data());
   AliAnalysisDataContainer *cTritonTree      = mgr->CreateContainer("reducedTree_Triton",      TTree::Class(),AliAnalysisManager::kOutputContainer,Filename.Data());
   AliAnalysisDataContainer *cHypertritonTree = mgr->CreateContainer("reducedTree_HyperTriton", TTree::Class(),AliAnalysisManager::kOutputContainer,Filename.Data());
   
   mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, cQA);
   mgr->ConnectOutput(task, 2, cOutput);
   mgr->ConnectOutput(task, 3, cHeliumTree);
   if (doTri)    mgr->ConnectOutput(task, 4, cTritonTree);
   if (doHypTri) mgr->ConnectOutput(task, 5, cHypertritonTree);
   
   
   AliAnalysisDataContainer *cTreeEventSelection = mgr->CreateContainer("TreeEventSelection", TTree::Class(),AliAnalysisManager::kOutputContainer,Filename.Data());
   mgr->ConnectOutput(task, 6, cTreeEventSelection);

   
   return task;
}
