AliAnalysisTaskStrangenessVsMultiplicity *AddTaskStrangenessVsMultiplicity( Bool_t lSaveV0 = kFALSE, Bool_t lSaveCascade = kTRUE, const TString lMasterJobSessionFlag = "")
{
// Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskStrangenessVsMultiplicity", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskStrangenessVsMultiplicity", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  // Create and configure the task
  AliAnalysisTaskStrangenessVsMultiplicity *taskAuxiliary = new AliAnalysisTaskStrangenessVsMultiplicity("taskAuxiliary");
  mgr->AddTask(taskAuxiliary);
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
   
  outputFileName += ":PWGLF_StrVsMult";
  if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
   
  Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );

  AliAnalysisDataContainer *coutputList = mgr->CreateContainer("cList",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );
  AliAnalysisDataContainer *coutputTree = mgr->CreateContainer("cTreeEvent",
  					       TTree::Class(),
						       AliAnalysisManager::kOutputContainer,
						       outputFileName );
  AliAnalysisDataContainer *coutputTreeV0 = mgr->CreateContainer("cTreeV0",
  					       TTree::Class(),
						       AliAnalysisManager::kOutputContainer,
						       outputFileName );
  AliAnalysisDataContainer *coutputTreeCascade = mgr->CreateContainer("cTreeCascade",
  					       TTree::Class(),
						       AliAnalysisManager::kOutputContainer,
						       outputFileName );
   
  //This one you should merge in file-resident ways...
  coutputTree->SetSpecialOutput();
  coutputTreeV0->SetSpecialOutput();
  coutputTreeCascade->SetSpecialOutput();

  //Recommendation: Tree as a single output slot
  mgr->ConnectInput (taskAuxiliary, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskAuxiliary, 1, coutputList);
  mgr->ConnectOutput(taskAuxiliary, 2, coutputTree);
  mgr->ConnectOutput(taskAuxiliary, 3, coutputTreeV0);
  mgr->ConnectOutput(taskAuxiliary, 4, coutputTreeCascade);
  
  return taskAuxiliary;
}   
