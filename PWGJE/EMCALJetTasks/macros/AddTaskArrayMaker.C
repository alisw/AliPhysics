/// \file AddTaskArrayMaker.C
/// \brief AddTask macro for the AliNanoAODArrayMaker class.
///
/// AddTask macro for the AliNanoAODArrayMaker class.
///
/// \author Markus Zimmermann
/// \date October 11, 2017

AliAnalysisTask * AddTaskArrayMaker(const char *output    = 0 , const char *outputPythia    = 0, const char *outputdata    = 0 ) {
  // Adds my task
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskArrayMaker", "No analysis manager to connect to.");
    return NULL;
  }  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskArrayMaker", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  // Configure analysis
  //===========================================================================

  AliNanoAODArrayMaker* task = new AliNanoAODArrayMaker("maker");
  task->SetOutputArrayName(output);
  task->SetOutputArrayPythiaName(outputPythia);
  task->SetOutputArrayDataName(outputdata);

  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  const char* containerName = "list";
  const char* folderName = "Array_Maker";
  TString fileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(containerName, TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", fileName.Data(), folderName));
  mgr->ConnectOutput(task, 1, coutput);

  return task;


}
