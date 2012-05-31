AliHMPIDPIDEfficiencyV0AnalysisTask *AddTaskHMPIDPIDEfficiencyV0(Bool_t useMC=kTRUE)
{
// Creates a HMPID task, configures it and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskHMPIDPIDEfficiencyV0", "No analysis manager to connect to.");
      return NULL;
   }

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskHMPIDPIDEfficiencyV0", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================

   AliHMPIDPIDEfficiencyV0AnalysisTask *hmpTask = new AliHMPIDPIDEfficiencyV0AnalysisTask("HMPIDPIDEfficiencyV0AnalysisTask");
   hmpTask->SetDebugLevel(0);
   hmpTask->SelectCollisionCandidates();
   hmpTask->SetUseMC(useMC);
   mgr->AddTask(hmpTask);

   AliAnalysisDataContainer *cout_hmpid= mgr->CreateContainer("HmpidOutput", TList::Class(),AliAnalysisManager::kOutputContainer,
                                           AliAnalysisManager::GetCommonFileName());
   AliAnalysisDataContainer *cout_tree = mgr->CreateContainer("HmpidTree", TTree::Class(),AliAnalysisManager::kOutputContainer,
                                           AliAnalysisManager::GetCommonFileName());

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (hmpTask, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (hmpTask, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (hmpTask, 1, cout_hmpid);
   mgr->ConnectOutput (hmpTask, 2, cout_tree);

   return hmpTask;
}
