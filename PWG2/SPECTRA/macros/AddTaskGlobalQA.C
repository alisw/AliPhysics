AliAnalysisTaskGlobalQA *AddTaskGlobalQA()
{
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskGlobalQA", "No analysis manager to connect to.");
      return NULL;
   }

   AliAnalysisTaskGlobalQA *taskGlobalQA = new AliAnalysisTaskGlobalQA();
   mgr->AddTask(taskGlobalQA);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   outname += ":PWG2GlobalQA";
   if (lCollidingSystems) outputFileName += "_AA";
   else outputFileName += "_PP";
   if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("GlobalQA",
							     TObjArray::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );

   mgr->ConnectInput (taskGlobalQA, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskGlobalQA, 1, coutput1);
   return taskGlobalQA;
}

