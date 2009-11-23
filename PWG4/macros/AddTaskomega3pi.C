AliAnalysisTaskOmegaPi0PiPi* AddTaskomega3pi()
{
  // Creates an omega(782) --> pi0 pi+ pi- analysis task, 
  // configures it and adds it to the analysis manager.

  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskomega3pi", "No analysis manager to connect to.");
    return NULL;
  }  

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskomega3pi", "This task requires an input event handler");
    return NULL;
  }
  
  // Add task
  AliAnalysisTaskOmegaPi0PiPi *omegaTask = new AliAnalysisTaskOmegaPi0PiPi("OmegaPi0PiPi");
  mgr->AddTask(omegaTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  //AliAnalysisDataContainer *coutput = mgr->CreateContainer("histos",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");


  TString outputfile = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("omega3pi", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:omega3pi",outputfile.Data()));
  
  // Connect input/output
  mgr->ConnectInput(omegaTask, 0, cinput);
  mgr->ConnectOutput(omegaTask, 1, coutput);
  
  return omegaTask;
}


