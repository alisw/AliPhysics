AliAnalysisTaskTOFTrigger *AddTaskTOFTrigger(const char *name,Float_t lowpt,Float_t highpt,Float_t highmult,TString trgcls,Int_t nBCs,Bool_t useEVS,Int_t cutSet, Float_t maxErr, Float_t mintof,Float_t maxtof, const char* suffix = ""){

  
  //--- get the current analysis manager ---//
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
      Error("AddTask_TOFTrigger", "No analysis manager found.");
      return 0;
   }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTask_TOFTrigger", "This task requires an input event handler");
    return 0;
  }
  
  TString combinedName;
  combinedName.Form("%s%s",name,suffix);
   	  
  // Create tasks
  AliAnalysisTaskTOFTrigger *task = new AliAnalysisTaskTOFTrigger(name,lowpt,highpt,highmult,trgcls,nBCs,useEVS,cutSet,maxErr,mintof,maxtof);
  mgr->AddTask(task);


   // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(combinedName, TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TOFTrig", AliAnalysisManager::GetCommonFileName())); 
  
  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

return task;
}
