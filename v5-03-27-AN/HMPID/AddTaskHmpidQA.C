AliHMPIDTaskQA *AddTaskHmpidQA(Bool_t useMC=kTRUE)
{
// Creates a HMPID task, configures it and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskHmpidQA", "No analysis manager to connect to.");
      return NULL;
   }

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskHmpidQA", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================

   AliHMPIDTaskQA *hmpTask = new AliHMPIDTaskQA("HMPIDTaskQA");
   hmpTask->SetDebugLevel(0);
   hmpTask->SelectCollisionCandidates();
   hmpTask->SetUseMC(useMC);
   mgr->AddTask(hmpTask);

   AliAnalysisDataContainer *cout_hmpid= mgr->CreateContainer("HmpidQA", TList::Class(),AliAnalysisManager::kOutputContainer,
                                           Form("%s:HmpidQA",AliAnalysisManager::GetCommonFileName()));
   AliAnalysisDataContainer *cout_tree = mgr->CreateContainer("HmpidTree", TTree::Class(),AliAnalysisManager::kOutputContainer,
                                           Form("%s:HmpidQA",AliAnalysisManager::GetCommonFileName()));

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (hmpTask, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (hmpTask, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (hmpTask, 1, cout_hmpid);
   mgr->ConnectOutput (hmpTask, 2, cout_tree);

   return hmpTask;
}