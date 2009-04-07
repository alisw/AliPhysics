AliAnalysisTaskParticleCorrelation *AddTaskPi0()
{
// Creates a Pi0->gg identification task and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTaskGammaHadronCorr", "No analysis manager to connect to.");
      return NULL;
   }   

   // This task requires an ESD input handler and an AOD output handler.
   // Check this using the analysis manager.
   //===============================================================================
   TString type = mgr->GetInputEventHandler()->GetDataType();
   if (!type.Contains("ESD")) {
      Error("AddTaskESDFilter", "ESD filtering task needs the manager to have an ESD input handler.");
      return NULL;
   }   
   // Check if AOD output handler exist.
   AliAODHandler *aod_h = (AliAODHandler*)mgr->GetOutputEventHandler();
   if (!aod_h) {
      Error("AddTaskESDFilter", "ESD filtering task needs the manager to have an AOD output handler.");
      return NULL;
   }   
   
   // Create the task, add it to the manager and configure it.
   //===========================================================================
   AliAnalysisTaskParticleCorrelation * taskpartcorr = new AliAnalysisTaskParticleCorrelation ("PWG4-Pi0");
   mgr->AddTask(taskpartcorr);
   taskpartcorr->SetConfigFileName("ConfigAnalysisPi0");
   // Output histograms list for particle correlation analysis                       
   AliAnalysisDataContainer *cout_partcorr = mgr->CreateContainer("pi0histos", TList::Class(),
			      AliAnalysisManager::kOutputContainer, "pi0histos.root");
   
   // Connect to data containers
   mgr->ConnectInput  (taskpartcorr,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (taskpartcorr,  1, cout_partcorr );
   return taskpartcorr;
}
   
