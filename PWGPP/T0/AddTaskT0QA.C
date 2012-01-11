AliAnalysisTaskSE* AddTaskT0QA()
{
  // Creates a QA task to check T0 data
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskT0QA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskT0QA", "This task requires an input event handler");
    return NULL;
  }
   TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
   // Configure analysis
   //===========================================================================
   
 
 
   AliT0AnalysisTaskQA* task = new AliT0AnalysisTaskQA("AliAnaT0QA");
   mgr->AddTask(task);
  
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("QAT0chists", TObjArray::Class(), 
							   AliAnalysisManager::kOutputContainer, Form("%s:T0_Performance",
												      mgr->GetCommonFileName()));

  
   mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (task, 1, coutput);

   return task;   
}


