AliAnalysisTaskCheckV0 *AddTaskCheckV0(Short_t lCollidingSystems=0  /*0 = pp, 1 = AA*/)
{
// Creates, configures and attaches to the train a V0 check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskCheckV0", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskCheckV0", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Create and configure the task
	AliAnalysisTaskCheckV0 *taskcheckv0 = new AliAnalysisTaskCheckV0("TaskCheckV0");
   taskcheckv0->SetCollidingSystems(lCollidingSystems);
   taskcheckv0->SetAnalysisType(type);
   mgr->AddTask(taskcheckv0);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   TString outname = "PP";
   if (lCollidingSystems) outname = "AA";
   if (mgr->GetMCtruthEventHandler()) outname += "-MC-";
   outname += "V0List.root";
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clistV0",
								   TList::Class(),
								   AliAnalysisManager::kOutputContainer,
								   outname );
                           
	mgr->ConnectInput(taskcheckv0, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskcheckv0, 1, coutput1);
   return taskcheckv0;
}   
