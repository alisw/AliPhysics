AliAnalysisTaskJetSpectrum *AddTaskJetSpectrum()
{
// Creates a jet fider task, configures it and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskJetSpectrum", "No analysis manager to connect to.");
      return NULL;
   }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskJetSpectrum", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================
   
   AliAnalysisTaskJetSpectrum* pwg4spec = new  AliAnalysisTaskJetSpectrum("Jet Spectrum");
      
   // or a config file
   pwg4spec->SetAnalysisType(AliAnalysisTaskJetSpectrum::kAnaMC);
   //      if(iAODanalysis)pwg4spec->SetAODInput(kTRUE);
   // pwg4spec->SetDebugLevel(11); 
   //      pwg4spec->SetBranchRec("jetsMC"); 
   //      pwg4spec->SetBranchGen("jetsMC"); 
   mgr->AddTask(pwg4spec);



      
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_Spec = mgr->CreateContainer("pwg4spec", TList::Class(),AliAnalysisManager::kOutputContainer,"pwg4spec.root");

   mgr->ConnectInput  (pwg4spec, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (pwg4spec, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (pwg4spec,  1, coutput1_Spec );
   
   return pwg4spec;
}

AliAnalysisTaskJetSpectrum *AddTaskJetSpectrum(AliAnalysisManager* mgr = 0,AliAnalysisDataContainer *cinput = 0)
{
  // This is only for running on PROOF with the old root version 5-22-00 
  // and the older version of the AF

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskJetSpectrum", "No analysis manager to connect to.");
      return NULL;
   }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskJetSpectrum", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================
   AliAnalysisTaskJetSpectrum* pwg4spec = new  AliAnalysisTaskJetSpectrum("Jet Spectrum");
   pwg4spec->SetAnalysisType(AliAnalysisTaskJetSpectrum::kAnaMC);
   //      if(iAODanalysis)pwg4spec->SetAODInput(kTRUE);
   //       pwg4spec->SetDebugLevel(11); 
   //      pwg4spec->SetBranchRec("jetsMC"); 
   //      pwg4spec->SetBranchGen("jetsMC"); 
   mgr->AddTask(pwg4spec);

   AliAnalysisDataContainer *coutput1_Spec = mgr->CreateContainer("pwg4spec", TList::Class(),AliAnalysisManager::kOutputContainer,"pwg4spec.root");

   // Dummy AOD output container for jet analysis (no client yet)
   c_aodSpec = mgr->CreateContainer("cAODjetSpec", TTree::Class(),
				    AliAnalysisManager::kExchangeContainer);
   mgr->ConnectInput  (pwg4spec,  0, cinput);    
   mgr->ConnectOutput (pwg4spec,  0, c_aodSpec );
   mgr->ConnectOutput (pwg4spec,  1, coutput1_Spec );

   return pwg4spec;

}  
