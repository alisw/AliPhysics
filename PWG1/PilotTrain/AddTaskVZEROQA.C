AliAnalysisTaskSE* AddTaskVZEROQA(Int_t runNumber)
{
  // Creates a QA task exploiting simple symmetries phi, eta +/-, charge ...
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskQAsym", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTasQAsym", "This task requires an input event handler");
    return NULL;
  }
   TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
   // Configure analysis
   //===========================================================================
   
 
 
   AliAnaVZEROQA* task = new AliAnaVZEROQA("AliAnaVZEROQA");
   mgr->AddTask(task);
  
   AliAnalysisDataContainer *cout  = mgr->CreateContainer("QAVZEROHists",TList::Class(),
							  AliAnalysisManager::kOutputContainer, Form("%s:VZERO_Performance", 
												     mgr->GetCommonFileName()));

   mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (task, 1, cout);

   return task;
   
  
}


