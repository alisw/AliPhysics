AliAnalysisTaskCaloFilter * AddTaskCaloFilter(){

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCaloFilter", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskCaloFilter", "This task requires an input event handler");
    return NULL;
  }
  
  AliAnalysisTaskCaloFilter * filter = new AliAnalysisTaskCaloFilter("CaloFilter");
  filter->SetCaloFilter(AliAnalysisTaskCaloFilter::kBoth); //kPHOS, kEMCAL or kBoth
  filter->SwitchOffClusterCorrection();
  
  filter->SetVzCut(10.);
  filter->SetEnergyCut(10.);
  filter->SetNcellsCut(2);
  
  filter->SwitchOnFillTracks();
  filter->SwitchOnFillAODFile();
  
  filter->PrintInfo(); 
  mgr->AddTask(filter);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();
    
  mgr->ConnectInput  (filter, 0, cinput1);
  mgr->ConnectOutput (filter, 0, coutput1 );
  
  return filter;

}