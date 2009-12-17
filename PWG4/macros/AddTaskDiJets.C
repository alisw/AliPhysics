AliAnalysisTaskDiJets *AddTaskDiJets()
{
// Creates a dijet task, configures it and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskJets", "No analysis manager to connect to.");
      return NULL;
   }

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskDiJets", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================

   AliAnalysisTaskDiJets *dijetana = new AliAnalysisTaskDiJets("DiJetAnalysis");
   dijetana->SetDebugLevel(10);
   mgr->AddTask(dijetana);

   AliAnalysisDataContainer *cout_dijet = mgr->CreateContainer("DiJet", TList::Class(),AliAnalysisManager::kOutputContainer,
     Form("%s:PWG4_DiJet",AliAnalysisManager::GetCommonFileName()));

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (dijetana, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (dijetana, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (dijetana, 1, cout_dijet);

   return dijetana;
}
