AliAnalysisTaskCaloConv * AddTaskCaloConv(){
  //Macro to add class CaloConv (conversion+calorimeters pi0 analysis) to train
  //Argument is the path to the PHOS recalibration parameters (file with OCDB entry)
  //Default path to the file with unit recalibration == no recalibnration
  //If file does not exist, no recalibration histograms will be filled

  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCaloConv", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskCaloConv", "This task requires an input event handler");
    return NULL;
  }

  // Add task
  AliAnalysisTaskCaloConv *task = new AliAnalysisTaskCaloConv("CaloConv");
  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  TString outputfile = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("CaloConv", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:PWG4_CaloConv",outputfile.Data()));
  mgr->ConnectOutput(task, 1, coutput);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("CFCaloConv", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:PWG4_CFCaloConv",outputfile.Data()));
  mgr->ConnectOutput(task, 2, coutput2);

  return task ;

}
