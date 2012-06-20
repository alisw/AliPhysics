AliAnalysisTaskJetsHMPID *AddTaskJetsHMPID(Char_t *jb="jets", Char_t *bb="")
{
// Creates a HMPID task for jet chemistry analysis, configures it and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskHMPID", "No analysis manager to connect to.");
      return NULL;
   }

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskJetsHMPID", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================

   AliAnalysisTaskJetsHMPID *hmpJetsTask = new AliAnalysisTaskJetsHMPID("JetsHMPIDAnalysisTask");
   hmpJetsTask->SetDebugLevel(0);
   hmpJetsTask->SelectCollisionCandidates();
   hmpJetsTask->SetJetBranch(jb);
   hmpJetsTask->SetBkgBranch(bb);
   hmpJetsTask->SetJetPtCut(10.);
   mgr->AddTask(hmpJetsTask);

   AliAnalysisDataContainer *cout_hmpid= mgr->CreateContainer("JetsHmpidOutput", TList::Class(),AliAnalysisManager::kOutputContainer,
                                           AliAnalysisManager::GetCommonFileName());
   AliAnalysisDataContainer *cout_tree = mgr->CreateContainer("JetsHmpidTree", TTree::Class(),AliAnalysisManager::kOutputContainer,
                                           AliAnalysisManager::GetCommonFileName());

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (hmpJetsTask, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (hmpJetsTask, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (hmpJetsTask, 1, cout_hmpid);
   mgr->ConnectOutput (hmpJetsTask, 2, cout_tree);

   return hmpJetsTask;
}
