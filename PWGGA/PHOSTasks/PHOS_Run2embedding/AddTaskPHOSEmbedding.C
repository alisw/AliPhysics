class AliAnalysisDataContainer;
AliPHOSEmbeddingRun2* AddTaskPHOSEmbedding(const char* name = "PHOSEmbedding")
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));

  AliPHOSEmbeddingRun2* task = new AliPHOSEmbeddingRun2(Form("%sTask", name));
  task->SelectCollisionCandidates(AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kINT7);
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  AliAODHandler* aodHandler = new AliAODHandler();
  aodHandler->SetOutputFileName("AliAODout.root");
  mgr->SetOutputEventHandler(aodHandler);
  AliAnalysisDataContainer* coutput =
    mgr->CreateContainer("output0", TTree::Class(), AliAnalysisManager::kOutputContainer);
  mgr->ConnectOutput(task, 0, coutput);

  AliPHOSEmbedggHBT* task1 = new AliPHOSEmbedggHBT(Form("%s_kCentral", name));
  task1->SelectCollisionCandidates(AliVEvent::kCentral);
  mgr->AddTask(task1);
  AliAnalysisDataContainer* coutput1 =
    mgr->CreateContainer(Form("%s_kCentral", name), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectInput(task1, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task1, 1, coutput1);

  AliPHOSEmbedggHBT* task2 = new AliPHOSEmbedggHBT(Form("%s_kSemiCentral", name));
  task2->SelectCollisionCandidates(AliVEvent::kSemiCentral);
  mgr->AddTask(task2);
  AliAnalysisDataContainer* coutput2 = mgr->CreateContainer(Form("%s_kSemiCentral", name), TList::Class(),
                                                            AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectInput(task2, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task2, 1, coutput2);

  AliPHOSEmbedggHBT* task3 = new AliPHOSEmbedggHBT(Form("%s_kINT7", name));
  task3->SelectCollisionCandidates(AliVEvent::kINT7);
  mgr->AddTask(task3);
  AliAnalysisDataContainer* coutput3 =
    mgr->CreateContainer(Form("%s_kINT7", name), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectInput(task3, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task3, 1, coutput3);

  return task;
}
