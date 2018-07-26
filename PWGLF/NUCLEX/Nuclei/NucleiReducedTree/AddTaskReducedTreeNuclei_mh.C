//______________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTask *AddTaskReducedTreeNuclei_mh()  {

   //Get the current analysis manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) 
   {
      Error("No analysis manager found.");
      return 0;
   }
   
   AliAnalysisTaskReducedTreeNuclei *task = new AliAnalysisTaskReducedTreeNuclei("TaskNuclei");
   task->SelectCollisionCandidates(AliVEvent::kINT7);
   task->useTritonCandidates(kTRUE);   
   mgr->AddTask(task);
   
   TString Filename = "AnalysisResults.root"; //mgr->GetCommonFileName();
   AliAnalysisDataContainer *cQA = mgr->CreateContainer("QAHistograms", TList::Class(),AliAnalysisManager::kOutputContainer,Filename.Data());
   Filename	+= ":Trees";
   AliAnalysisDataContainer *cTreeEventSelection = mgr->CreateContainer("TreeEventSelection", TTree::Class(),AliAnalysisManager::kOutputContainer,Filename.Data());
   AliAnalysisDataContainer *cHeliumTree = mgr->CreateContainer("reducedTree_Helium", TTree::Class(),AliAnalysisManager::kOutputContainer,Filename.Data());
   AliAnalysisDataContainer *cHypertritonTree = mgr->CreateContainer("reducedTree_HyperTriton", TTree::Class(),AliAnalysisManager::kOutputContainer,Filename.Data());
   mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());

   mgr->ConnectOutput(task, 1, cQA);
   mgr->ConnectOutput(task, 2, cTreeEventSelection);
   mgr->ConnectOutput(task, 3, cHeliumTree);
   mgr->ConnectOutput(task, 4, cHypertritonTree);

   return task;
}
