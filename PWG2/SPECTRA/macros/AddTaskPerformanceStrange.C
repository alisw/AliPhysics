AliAnalysisTaskPerformanceStrange *AddTaskPerformanceStrange(Short_t lCollidingSystems=0,  /*0 = pp, 1 = AA*/
				       TString lAnalysisCut="no" )
{
// Creates, configures and attaches to the train a strangeness task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskPerformanceStrange", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskPerformanceStrange", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Create and configure the task
	AliAnalysisTaskPerformanceStrange *taskperformancestrange = new AliAnalysisTaskPerformanceStrange("TaskPerformanceStrange");
   taskperformancestrange->SetCollidingSystems(lCollidingSystems);
   taskperformancestrange->SetAnalysisType(type);
   taskperformancestrange->SetAnalysisCut(lAnalysisCut);
   mgr->AddTask(taskperformancestrange);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   outputFileName += ":PWG2Strange";
   if (lCollidingSystems) outputFileName += "_AA";
   else outputFileName += "_PP";
   if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clistPerformanceStrange",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );
                           
   mgr->ConnectInput(taskperformancestrange, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskperformancestrange, 1, coutput1);
   return taskperformancestrange;
}   
