void AddTaskFragmentationTriggered()
//AliAnalysisTaskFragmentationTriggered *AddTaskFragmentationTriggered()
{
  const char *name = "fragmentationtriggered";

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  //Min Bias sample
  AliAnalysisTaskFragmentationTriggered *task = new AliAnalysisTaskFragmentationTriggered("fragmentationtriggered");
  task->SetTrigger(0);
  mgr->AddTask(task);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->GetCommonOutputContainer();

  AliAnalysisDataContainer *clist =
    mgr->CreateContainer(Form("list_%s", name), TList::Class(), AliAnalysisManager::kOutputContainer,
                         Form("%s:fragmentationtriggered", AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectInput(task, 0, cinput);

  if (coutput)
    mgr->ConnectOutput(task, 0, coutput);

  mgr->ConnectOutput(task, 1, clist);


  //TRD triggered sample

 AliAnalysisTaskFragmentationTriggered *task_trg = new AliAnalysisTaskFragmentationTriggered("fragmentationtriggered_trg");

  task_trg->SetTrigger(1);
  mgr->AddTask(task_trg);


  name = "fragmentationtriggered_trg";
  AliAnalysisDataContainer *clist_trg =
    mgr->CreateContainer(Form("list_%s", name), TList::Class(), AliAnalysisManager::kOutputContainer,
                         Form("%s:fragmentationtriggered", AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectInput(task_trg, 0, cinput);

  if (coutput)
    mgr->ConnectOutput(task_trg, 0, coutput);

  mgr->ConnectOutput(task_trg, 1, clist_trg);

  return;
}
