AliAnalysisTaskJetSpectrum2 *AddTaskJetSpectrum2()
{
// Creates a jet fider task, configures it and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskJetSpectrum2", "No analysis manager to connect to.");
      return NULL;
   }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskJetSpectrum2", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================
   
   AliAnalysisTaskJetSpectrum2* pwg4spec = new  AliAnalysisTaskJetSpectrum2("Jet Spectrum");
      
   // or a config file
   pwg4spec->SetAnalysisType(AliAnalysisTaskJetSpectrum2::kAnaMC);
   //      if(iAODanalysis)pwg4spec->SetAODInput(kTRUE);
   // pwg4spec->SetDebugLevel(11); 
   pwg4spec->SetBranchGen("jetsMC"); 
   pwg4spec->SetBranchRec("jetsAOD"); 
   mgr->AddTask(pwg4spec);



      
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_Spec = mgr->CreateContainer("pwg4spec2", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_spec2",AliAnalysisManager::GetCommonFileName()));

   mgr->ConnectInput  (pwg4spec, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (pwg4spec, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (pwg4spec,  1, coutput1_Spec );
   
   return pwg4spec;
}
