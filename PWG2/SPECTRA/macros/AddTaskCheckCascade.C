AliAnalysisTaskCheckCascade *AddTaskCheckCascade(Short_t lCollidingSystems=0  /*0 = pp, 1 = AA*/)
{
// Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskCheckCascade", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskCheckCascade", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Create and configure the task
	AliAnalysisTaskCheckCascade *taskcheckcascade = new AliAnalysisTaskCheckCascade("TaskCheckCascade");
   taskcheckcascade->SetCollidingSystems(lCollidingSystems);
   taskcheckcascade->SetAnalysisType(type);
   mgr->AddTask(taskcheckcascade);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   TString outname = "PP";
   if (lCollidingSystems) outname = "AA";
   if (mgr->GetMCtruthEventHandler()) outname += "-MC-";
   outname += "CascadeList.root";
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clistCasc",
								   TList::Class(),
								   AliAnalysisManager::kOutputContainer,
								   outname );
                           
	mgr->ConnectInput(taskcheckcascade, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskcheckcascade, 1, coutput1);
   return taskcheckcascade;
}   
