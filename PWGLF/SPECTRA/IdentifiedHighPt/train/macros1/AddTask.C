/*
  Last update (new version): 
  20 aug 2015: clean up
  22 aug 2015: match pre defined output

*/

AliAnalysisTask* AddTask(const Char_t* taskname)
{
  
  // Creates a pid task and adds it to the analysis manager
  
  // Get the pointer to the existing analysis manager via the static
  //access methodh
  //=========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHighPtDeDx", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // Check the analysis type using the event handlers connected to the
  // analysis manager The availability of MC handler can also be
  // checked here.
  // =========================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskHighPtDeDx", "This task requires an input event handler");
    return NULL;
  }  
  
  // Create task  
  AliAnaTaskV0EffDecomposition* taskV0EffDecomposition 
    = new AliAnaTaskV0EffDecomposition(taskname);

  mgr->AddTask(taskV0EffDecomposition);
  
  // Create ONLY the output containers for the data produced by the
  // task.  Get and connect other common input/output containers via
  // the manager as below
  //=======================================================================
  

  Char_t outFileName[256]={0};
  sprintf(outFileName, "AnalysisResults.root");

  
  AliAnalysisDataContainer* cout_taskV0EffDecomposition = 
    mgr->CreateContainer("outputV0EffDecomposition", TList::Class(), AliAnalysisManager::kOutputContainer, outFileName);
  mgr->ConnectInput (taskV0EffDecomposition, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskV0EffDecomposition, 1, cout_taskV0EffDecomposition);

    
  // Return task pointer at the end
  return taskV0EffDecomposition;

}

