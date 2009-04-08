AliAnalysisTaskStrange *AddTaskStrange(Short_t lCollidingSystems=0,  /*0 = pp, 1 = AA*/
                                       const char *optCuts="")
{
// Creates, configures and attaches to the train a strangeness task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskStrange", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskStrange", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Create and configure the task
	AliAnalysisTaskStrange *taskstrange = new AliAnalysisTaskStrange("TaskStrange", optCuts);
   taskstrange->SetCollidingSystems(lCollidingSystems);
   taskstrange->SetAnalysisType(type);
   mgr->AddTask(taskstrange);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   TString outname = "PP";
   if (lCollidingSystems) outname = "AA";
   if (mgr->GetMCtruthEventHandler()) outname += "-MC-";
   outname += "StrangeList.root";
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clistStrange",
								   TList::Class(),
								   AliAnalysisManager::kOutputContainer,
								   outname );
                           
	mgr->ConnectInput(taskstrange, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskstrange, 1, coutput1);
   return taskstrange;
}   
