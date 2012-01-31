AliVZEROEPSelectionTask *AddTaskVZEROEPSelection()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskVZEROEPSelection", "No analysis manager to connect to.");
    return NULL;
  }    
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskVZEROEPSelection", "This task requires an input event handler");
    return NULL;
  }

  AliVZEROEPSelectionTask *task = new AliVZEROEPSelectionTask("AliVZEROEPSelectionTask");
  //  task->UseVZEROCentrality(); // Optional line to swith from SPD to VZERO centrality estimator
  mgr->AddTask(task);
  
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  return task;
}   
