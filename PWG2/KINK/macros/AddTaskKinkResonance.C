AliResonanceKinkPID *AddTaskKinkResonance(Short_t lCollidingSystems=0  /*0 = pp, 1 = AA*/)
{
// Creates, configures and attaches to the train a V0 check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskKinkResonance", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskKinkResonance", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
   if (type != "ESD") {
      ::Error("AddTaskKinkResonance", "This task needs an ESD input handler");
      return NULL;
   }   
   if (!mgr->GetMCtruthEventHandler()) {
      ::Error("AddTaskKinkResonance", "This task needs an MC handler");
      return NULL;
   }

   // Create and configure the task
	AliResonanceKinkPID *taskkinkresonance = new AliResonanceKinkPID("TaskResKinkPID");
   mgr->AddTask(taskkinkresonance);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   TString outname = "PP";
   if (lCollidingSystems) outname = "AA";
   if (mgr->GetMCtruthEventHandler()) outname += "-MC-";
   outname += "KinkResonanceList.root";
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("KinkResPID",
								   TList::Class(),
								   AliAnalysisManager::kOutputContainer,
								   outname );
                           
	mgr->ConnectInput(taskkinkresonance, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskkinkresonance, 1, coutput1);
   return taskkinkresonance;
}   
