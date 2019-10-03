AliAnalysisTaskExtractCascade *AddTaskExtractCascade( Bool_t lSwitchIsNuclear     = kFALSE, 
                                                                  Bool_t lSwitchINT7 = kFALSE,
                                                                  const TString lMasterJobSessionFlag = "")
{
// Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskExtractCascade", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskExtractCascade", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Create and configure the task
	 AliAnalysisTaskExtractCascade *taskextract = new AliAnalysisTaskExtractCascade("taskextract");

   //Configuration
   taskextract -> SetIsNuclear     ( lSwitchIsNuclear     );
   taskextract -> SetINT7Trigger   ( lSwitchINT7          );

   mgr->AddTask(taskextract);

   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   
   outputFileName += ":PWGLFExtractCascade";
   //if (lSwitchIsNuclear) outputFileName += "_AA";
   outputFileName += "_PP";
   if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
   //if(lMasterJobSessionFlag.Length()) outputFileName += lMasterJobSessionFlag.Data();
   
   Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );

   AliAnalysisDataContainer *coutputList = mgr->CreateContainer("clist",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );
   AliAnalysisDataContainer *coutputTreeCascade = mgr->CreateContainer("cTreeCascade",
							     TTree::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );
   
   //This one you should merge in file-resident ways...
   coutputTreeCascade->SetSpecialOutput();

   //Recommendation: Tree as a single output slot
   mgr->ConnectInput( taskextract, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskextract, 1, coutputList);
   mgr->ConnectOutput(taskextract, 2, coutputTreeCascade);
   
   return taskextract;
}   
