AliAnalysisTaskExtractPerformanceV0 *AddTaskExtractPerformanceV0( Bool_t lSwitchIsNuclear     = kTRUE, 
                                                                  Bool_t lSwitchIsLowEnergyPP = kFALSE,
                                                                  Bool_t lSwitchUseOnTheFly   = kFALSE, 
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
   taskv0extractperformance -> SetIsLowEnergyPP ( lSwitchIsLowEnergyPP );
   taskv0extractperformance -> SetUseOnTheFly   ( lSwitchUseOnTheFly   );

   mgr->AddTask(taskv0extractperformance);

   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   
   outputFileName += ":PWG2CheckPerformanceLambda";

   else outputFileName += "_PP";
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

   mgr->ConnectInput( taskv0extractperformance, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskv0extractperformance, 1, coutputList);
   mgr->ConnectOutput(taskv0extractperformance, 2, coutputTree);
   
   return taskv0extractperformance;
}   
