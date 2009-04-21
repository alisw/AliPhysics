AliAnalysisTaskSEVertexingHF *AddTaskVertexingHF() {
  //
  // Creates a task for heavy flavour vertexing and adds it to the analysis manager.
  // andrea.dainese@lnl.infn.it
  //

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskVertexingHF", "No analysis manager to connect to.");
    return NULL;
  }   
   
  // This task requires an ESD or AOD input handler and an AOD output handler.
  // Check this using the analysis manager.
  //===============================================================================
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddTaskVertexingHF", "HF vertexing task needs the manager to have an ESD or AOD input handler.");
    return NULL;
  }   
  // Check if AOD output handler exist.
  AliAODHandler *aodh = (AliAODHandler*)mgr->GetOutputEventHandler();
  if (!aodh) {
    ::Error("AddTaskVertexingHF", "HF vertexing task needs the manager to have an AOD output handler.");
    return NULL;
  }   
  
  // Create the task, add it to the manager and configure it.
  //===========================================================================
  AliAnalysisTaskSEVertexingHF *hfTask = new AliAnalysisTaskSEVertexingHF("vertexing HF");
  mgr->AddTask(hfTask);

  //
  // Create containers for input/output
  mgr->ConnectInput(hfTask,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(hfTask,0,mgr->GetCommonOutputContainer());

  return hfTask;
}
