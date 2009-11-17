AliAnalysisTaskCheckPerformanceCascade *AddTaskCheckPerformanceCascade(Short_t lCollidingSystems=0  /*0 = pp, 1 = AA*/ )
{
// Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskCheckPerformanceCascade", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskCheckPerformanceCascade", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Create and configure the task
   AliAnalysisTaskCheckPerformanceCascade *taskCheckPerfCascade = new AliAnalysisTaskCheckPerformanceCascade("TaskCheckPerfCascade");
   taskCheckPerfCascade->SetCollidingSystems(lCollidingSystems);
   taskCheckPerfCascade->SetAnalysisType(type);
   mgr->AddTask(taskCheckPerfCascade);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   outname += ":PWG2CheckPerformanceCascade";
   if (lCollidingSystems) outputFileName += "_AA";
   else outputFileName += "_PP";
   if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clistCascMC",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );

   mgr->ConnectInput( taskCheckPerfCascade, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskCheckPerfCascade, 1, coutput1);
   return taskCheckPerfCascade;
}   
