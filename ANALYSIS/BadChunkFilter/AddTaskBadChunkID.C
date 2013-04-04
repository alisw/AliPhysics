AliAnalysisTaskBadChunkID *AddTaskBadChunkID( const TString lMasterJobSessionFlag = "")
{
// Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskBadChunkID", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskBadChunkID", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Create and configure the task
	 AliAnalysisTaskBadChunkID *taskbadchunkID = new AliAnalysisTaskBadChunkID("taskbadchunkID");
   mgr->AddTask(taskbadchunkID);

   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   
   outputFileName += ":BADCHUNKCHECK";
   //if (lSwitchIsNuclear) outputFileName += "_AA";
   //outputFileName += "_PP";
   if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
   //if(lMasterJobSessionFlag.Length()) outputFileName += lMasterJobSessionFlag.Data();
   
   Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );

   AliAnalysisDataContainer *coutputList = mgr->CreateContainer("clist",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );
   AliAnalysisDataContainer *coutputTree = mgr->CreateContainer("cTree",
							     TTree::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );

   //This one you should merge in file-resident ways...
   coutputTree->SetSpecialOutput();

   //Recommendation: Tree as a single output slot
   mgr->ConnectInput( taskbadchunkID, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskbadchunkID, 1, coutputList);
   mgr->ConnectOutput(taskbadchunkID, 2, coutputTree);
   
   return taskbadchunkID;
}   
