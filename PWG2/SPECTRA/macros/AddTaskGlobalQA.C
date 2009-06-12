AliAnalysisTaskGlobalQA *AddTaskGlobalQA()
{
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskGlobalQA", "No analysis manager to connect to.");
      return NULL;
   }

   AliAnalysisTaskGlobalQA *taskGlobalQA = new AliAnalysisTaskGlobalQA();
   mgr->AddTask(taskGlobalQA);

   AliAnalysisDataContainer *coutput1 = 
      mgr->CreateContainer("GlobalQA", TObjArray::Class(),
      AliAnalysisManager::kOutputContainer, "GlobalQA.root");

   mgr->ConnectInput (taskGlobalQA, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskGlobalQA, 1, coutput1);
   return taskGlobalQA;
}

