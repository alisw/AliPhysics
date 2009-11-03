AliAnalysisTaskMEVertexingHF *AddTaskHFMixing() {
  //
  // Creates a task for event mixing and adds it to the analysis manager.
  // r.romita@gsi.de
  //

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHFMixing", "No analysis manager to connect to.");
    return NULL;
  }   
   
  // This task requires  AOD input handler and an AOD output handler.
  // Check this using the analysis manager.
  //===============================================================================

  // Check if AOD output handler exist.
  AliAODHandler *aodh = (AliAODHandler*)mgr->GetOutputEventHandler();
  if (!aodh) {
    ::Error("AddTaskMixing", "HF vertexing task needs the manager to have an AOD output handler.");
    return NULL;
  }   
  
  // Create the task, add it to the manager and configure it.
  //===========================================================================
  AliAnalysisTaskMEVertexingHF *hfTask = new AliAnalysisTaskMEVertexingHF("mixing vertexing HF");
  mgr->AddTask(hfTask);

  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain",TChain::Class(), AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(hfTask,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(hfTask,0,mgr->GetCommonOutputContainer());

  return hfTask;
}
