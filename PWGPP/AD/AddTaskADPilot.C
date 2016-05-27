AliAnalysisTaskADPilot *AddTaskADPilot(){

  
  //--- get the current analysis manager ---//
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
      Error("AddTask_ADQA", "No analysis manager found.");
      return 0;
   }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTask_ADQA", "This task requires an input event handler");
    return 0;
  }
	
    
  // Create tasks
  AliAnalysisTaskADPilot *task = new AliAnalysisTaskADPilot("ESD");
  mgr->AddTask(task);

   // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("ADQAListHist", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ADQA", AliAnalysisManager::GetCommonFileName()));
 
  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

return task;
}
