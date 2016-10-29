AliEmcalCorrectionTask* AddTaskEmcalCorrectionTask(const char * suffix = "")
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalCorrectionTask", "No analysis manager to connect to.");
    return 0;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalCorrectionTask", "This task requires an input event handler");
    return 0;
  }
  
  TString name = "AliEmcalCorrectionTask";
  if (suffix != "") {
    name += TString::Format("_%s", suffix);
  }
  
  AliEmcalCorrectionTask* mgrTask = static_cast<AliEmcalCorrectionTask *>(mgr->GetTask(name.Data()));
  if (mgrTask) return mgrTask;
  
  // Create the task that manages the corrections
  AliEmcalCorrectionTask* correctionTask = new AliEmcalCorrectionTask(name.Data());

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(correctionTask);

  // Create containers for input/output
  AliAnalysisDataContainer* cInput = mgr->GetCommonInputContainer();
  
  TString outputContainerName(name);
  outputContainerName += "_histos";
  
  AliAnalysisDataContainer * cOutput = mgr->CreateContainer(outputContainerName.Data(),
                               TList::Class(),
                               AliAnalysisManager::kOutputContainer,
                               Form("%s", AliAnalysisManager::GetCommonFileName()));
  
  mgr->ConnectInput(correctionTask, 0, cInput);
  mgr->ConnectOutput(correctionTask, 1, cOutput);
  
  //TObjArray* cnt = mgr->GetContainers();

  return correctionTask;
}
