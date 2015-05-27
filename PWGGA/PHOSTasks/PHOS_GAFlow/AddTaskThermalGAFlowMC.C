AliAnalysisTaskThermalGAFlowMC* AddTaskThermalGAFlowMC(
   const char *outfilename    = "AnalysisOutput.root",
   const char *tag            = ""
)
{   
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskPHOSThermalGAFlowMC", "No analysis manager to connect to.");
    return NULL;
  }  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskPHOSThermalGAFlowMC", "This task requires an input event handler");
    return NULL;
  }
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString name(Form("ThermalGA_%s", tag));
  AliAnalysisTaskThermalGAFlowMC *PHOSGAtask = new AliAnalysisTaskThermalGAFlowMC(name);
  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(PHOSGAtask);
  // Create containers for input/output
  mgr->ConnectInput (PHOSGAtask, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *coGAFlow = mgr->CreateContainer(name,
    TList::Class(),
    AliAnalysisManager::kOutputContainer,
    outfilename);

  mgr->ConnectOutput(PHOSGAtask,1,coGAFlow);
  return PHOSGAtask;
}

//Note: Order of global variables matter!  They must match header, addtask, and cxx.
