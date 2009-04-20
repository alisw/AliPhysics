AliAnalysisTaskDiJets *AddTaskDiJets()
{
// Creates a jet fider task, configures it and adds it to the analysis manager.

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
   
   mgr->AddTask(dijetana);
      
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (dijetana, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (dijetana, 0, mgr->GetCommonOutputContainer());
   
   return dijetana;
}

AliAnalysisTaskDiJets *AddTaskDiJets(AliAnalysisManager* mgr,AliAnalysisDataContainer *cinput)
{
  // This is only for running on PROOF with the old root version 5-22-00 
  // and the older version of the AF

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
   mgr->AddTask(dijetana);

   //
   // Create containers for input/output
   AliAnalysisDataContainer *c_aod_dijet = mgr->CreateContainer("cAODdijet", TTree::Class(),AliAnalysisManager::kExchangeContainer);
   mgr->ConnectInput  (dijetana,  0, cinput  );
   mgr->ConnectOutput (dijetana,  0, c_aod_dijet);

   return dijetana;

}  
