AliAnalysisTaskPHOSPi0CalibSelection *AddTaskPHOSPi0Calib()
{
  // Creates a PartCorr task, configures it and adds it to the analysis manager.
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPartCorr", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPartCorr", "This task requires an input event handler");
    return NULL;
  }

   //TString dataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"  
   // Configure analysis
   //===========================================================================
    
   // Create task
   //===========================================================================

  AliAnalysisTaskPHOSPi0CalibSelection * pi0calib = new AliAnalysisTaskPHOSPi0CalibSelection ("PHOSPi0Calibration");
  //pi0calib->CopyAOD(kTRUE);
  mgr->AddTask(pi0calib);
  
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  
 AliAnalysisDataContainer *coutput = 
    mgr->CreateContainer("PHOSPi0Calib", TList::Class(), AliAnalysisManager::kOutputContainer, "PHOSPi0Calib.root");

  mgr->ConnectInput  (pi0calib, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (pi0calib, 1, coutput);
  
  return pi0calib;
}


