AliAnalysisTaskSE* AddTaskZDCQApA()
{
  // Creates a QA task to check ZDC data

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskZDCQA", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskZDCQA", "This task requires an input event handler");
    return NULL;
  }
   TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Configure analysis
   //===========================================================================
   //AliAnalysisTaskZDCpp* task = new AliAnalysisTaskZDCpp("AliAnaZDCQA");
   //AliAnalysisTaskZDCPbPb* task = new AliAnalysisTaskZDCPbPb("AliAnaZDCQA");
   AliAnalysisTaskZDCpA* task = new AliAnalysisTaskZDCpA("AliAnaZDCQA");
   mgr->AddTask(task);

   AliAnalysisDataContainer *cout  = mgr->CreateContainer("QAZDCHists",TList::Class(),
							  AliAnalysisManager::kOutputContainer, Form("%s:ZDC_Performance",
												     mgr->GetCommonFileName()));

   mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (task, 1, cout);

   return task;
}
