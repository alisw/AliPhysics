AliAnalysisTaskFragmentationFunction *AddTaskFragmentationFunction(UInt_t filterMask = 0,Int_t iPhysicsSelection)
{
// Creates a jet fider task, configures it and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskFragmentationFunction", "No analysis manager to connect to.");
      return NULL;
   }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskMyDyJets", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================
   
   AliAnalysisTaskFragmentationFunction* pwg4dijets = new  AliAnalysisTaskFragmentationFunction("Fragmentation Function Study");
      
   pwg4dijets->SetBranchGen("jetsMC"); 
   pwg4dijets->SetBranchRec("jets"); 
   //   pwg4dijets->SetBranchRec("jetsUA1AOD");
   pwg4dijets->SetLimitGenJetEta(0);
   pwg4dijets->SetFilterMask(filterMask); 
    if(iPhysicsSelection) pwg4dijets->SelectCollisionCandidates();
   mgr->AddTask(pwg4dijets);
   



      
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_FF = mgr->CreateContainer("PWG4_FF", TList::Class(),AliAnalysisManager::kOutputContainer,"PWG4_FF.root");

   mgr->ConnectInput  (pwg4dijets, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (pwg4dijets, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (pwg4dijets,  1, coutput1_FF );
   
   return pwg4dijets;
}
