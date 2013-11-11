AliAnalysisTaskQAV0AOD *AddTaskQAV0AOD( TString lCustomName = "QA", const TString lMasterJobSessionFlag = "")
{
// Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskQAV0AOD", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskQAV0AOD", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Create and configure the task
	 AliAnalysisTaskQAV0AOD *taskv0extract = new AliAnalysisTaskQAV0AOD("taskv0extract");

   mgr->AddTask(taskv0extract);

   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   
   outputFileName += ":PWGLFQAV0";
   outputFileName += "_";
   outputFileName += lCustomName.Data(); 
   
   Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );

   AliAnalysisDataContainer *coutputList = mgr->CreateContainer("clist",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );

   mgr->ConnectInput( taskv0extract, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskv0extract, 1, coutputList);
   
   return taskv0extract;
}   
