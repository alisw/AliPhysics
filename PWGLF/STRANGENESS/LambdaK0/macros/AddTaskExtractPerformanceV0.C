AliAnalysisTaskExtractPerformanceV0 *AddTaskExtractPerformanceV0( Bool_t lSwitchIsNuclear     = kFALSE, 
                                                                  Bool_t lSwitchINT7          = kFALSE,
                                                                  Bool_t lSwitchUseOnTheFly   = kFALSE,
                                                                  Bool_t lSwitchTakeAllTracks = kFALSE,  
                                                                  const TString lMasterJobSessionFlag = "")
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
	 AliAnalysisTaskExtractPerformanceV0 *taskv0extractperformance = new AliAnalysisTaskExtractPerformanceV0("taskv0extractperformance");

   //Configuration
   taskv0extractperformance -> SetIsNuclear     ( lSwitchIsNuclear     );
   taskv0extractperformance -> SetINT7Trigger   ( lSwitchINT7          );
   taskv0extractperformance -> SetUseOnTheFly   ( lSwitchUseOnTheFly   );
   taskv0extractperformance -> SetTakeAllTracks ( lSwitchTakeAllTracks );

   mgr->AddTask(taskv0extractperformance);

   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   
   outputFileName += ":PWGLFExtractPerformanceV0";
   //if (lSwitchIsNuclear) outputFileName += "_AA";
   outputFileName += "_PP";
   if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
   //if(lMasterJobSessionFlag.Length()) outputFileName += lMasterJobSessionFlag.Data();
   
   Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );

   AliAnalysisDataContainer *coutputList = mgr->CreateContainer("clistV0MC",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );
   AliAnalysisDataContainer *coutputTree = mgr->CreateContainer("cTreeMC",
							     TTree::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );
   
   //This one you should merge in file-resident ways...
   coutputTree->SetSpecialOutput();

   //Recommendation: Tree as a single output slot
   mgr->ConnectInput( taskv0extractperformance, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskv0extractperformance, 1, coutputList);
   mgr->ConnectOutput(taskv0extractperformance, 2, coutputTree);
   
   return taskv0extractperformance;
}   
