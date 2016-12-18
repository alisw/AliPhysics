AliAnalysisTaskExtractV0AODRun2* AddTaskExtractV0AODRun2( Bool_t lSwitchLowE     = kFALSE,
                                                   Bool_t lSwitchSaveAllInvMasses = kFALSE)
{
// Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTask2ExtractV0AOD", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTask2ExtractV0AOD", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
   // Create and configure the task
    AliAnalysisTaskExtractV0AODRun2* taskv0extract = new AliAnalysisTaskExtractV0AODRun2("V0Extract_PbPb_AOD");
    
   //Configuration
   taskv0extract -> SetIsNuclear     ( kTRUE                );
   taskv0extract -> SetIsLowEnergyPP ( kFALSE               );
   taskv0extract -> SetUseOnTheFly   ( kFALSE               );
    if(!lSwitchLowE) {
        taskv0extract -> SetTriggerMask   ( "kINT7"         );
    } else {
        taskv0extract -> SetTriggerMask   ( "kMB"         );
    }
    taskv0extract -> SetLowE         ( lSwitchLowE                );
    taskv0extract -> SetSaveAllInvMasses     ( lSwitchSaveAllInvMasses        );
    taskv0extract -> SetPreSelect     ( kTRUE        );

   mgr->AddTask(taskv0extract);

   TString outputFileName = AliAnalysisManager::GetCommonFileName();
        if (mgr->GetMCtruthEventHandler()) outputFileName.ReplaceAll(".root","_MC.root");
   
//   outputFileName += ":PWG2CheckLambda";
   //if (lCollidingSystems) outputFileName += "_AA_";
   //outputFileName += "_PP";
//   if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
   //if(lMasterJobSessionFlag.Length()) outputFileName += lMasterJobSessionFlag.Data();
   
    //Adjustments: two output files not allowed in the LEGO train framework at this moment
    Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );

   AliAnalysisDataContainer *coutputList = mgr->CreateContainer("clistV0",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );
   AliAnalysisDataContainer *coutputTree = mgr->CreateContainer("cTree",
							     TTree::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );
   
   //This one you should merge in file-resident ways...
   coutputTree->SetSpecialOutput();

   mgr->ConnectInput( taskv0extract, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskv0extract, 1, coutputList);
   mgr->ConnectOutput(taskv0extract, 2, coutputTree);
        
   return taskv0extract;
}   
