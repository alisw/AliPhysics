AliAnalysisTaskDiJets *AddTaskDiJets(Char_t *jb="jets")
{
// Creates a dijet task, configures it and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskDiJets", "No analysis manager to connect to.");
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

   AliAnalysisTaskDiJets *dijetana = new AliAnalysisTaskDiJets(Form("DiJetAnalysis_%s",jb));
   dijetana->SetDebugLevel(0);
//   dijetana->SetFillAOD(kTRUE);
   dijetana->SetJetBranch(jb);
   mgr->AddTask(dijetana);
   
   TString jbOut(jb);
   jbOut = jbOut(4,jbOut.Sizeof());
   jbOut.ToLower();

   AliAnalysisDataContainer *cout_dijet = mgr->CreateContainer(Form("dijets_%s",jbOut.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,
     Form("%s:PWG4_DiJets_%s",AliAnalysisManager::GetCommonFileName(),jbOut.Data()));

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (dijetana, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (dijetana, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (dijetana, 1, cout_dijet);

   return dijetana;
}
