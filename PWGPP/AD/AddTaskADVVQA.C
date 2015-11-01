AliAnalysisTaskADVVQA *AddTaskADVVQA(){

  
  //--- get the current analysis manager ---//
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
      Error("AddTask_ADVVQA", "No analysis manager found.");
      return 0;
   }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTask_ADVVQA", "This task requires an input event handler");
    return 0;
  }
	
    
  // Create tasks
  AliAnalysisTaskADVVQA *task = new AliAnalysisTaskADVVQA("ESD");
  mgr->AddTask(task);

   // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("ADVVQAListHist", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ADVVQA", AliAnalysisManager::GetCommonFileName()));  

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

return task;
}
