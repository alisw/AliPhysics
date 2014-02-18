AliAnalysisTaskVZEROEqFactorTask *AddTaskVZEROEqFactorTask( const TString lMasterJobSessionFlag = "")
{
// Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskExtractPerformanceV0", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskExtractPerformanceV0", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  // Create and configure the task
  AliAnalysisTaskVZEROEqFactorTask *taskVZEROAuxiliary = new AliAnalysisTaskVZEROEqFactorTask("taskVZEROAuxiliary");

  mgr->AddTask(taskVZEROAuxiliary);
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
   
  outputFileName += ":PWGLF_VZEROEqFactorTask";
  if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
   
  Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );

  AliAnalysisDataContainer *coutputList = mgr->CreateContainer("cList_VZERO",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );

  //Recommendation: Tree as a single output slot
  mgr->ConnectInput( taskVZEROAuxiliary, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskVZEROAuxiliary, 1, coutputList);
  
  return taskVZEROAuxiliary;
}   
