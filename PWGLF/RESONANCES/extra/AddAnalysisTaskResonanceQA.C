Bool_t AddAnalysisTaskResonanceQA()
{
   // retrieve analysis manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliAnalysisTaskResonanceQA *task = new AliAnalysisTaskResonanceQA("test");
   
   // initialize containers
   AliAnalysisDataContainer *out = mgr->CreateContainer("out", TList::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());

   // connect containers to AnalysisManager
   mgr->AddTask(task);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, out);
   
   return kTRUE;
}
