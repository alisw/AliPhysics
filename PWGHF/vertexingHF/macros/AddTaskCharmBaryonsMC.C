AliAnalysisTaskCharmBaryonsMC *AddTaskCharmBaryonsMC()
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCharmBaryonMC", "No analysis manager to connect to.");
    return NULL;
  }  


  //CREATE THE TASK

  printf("CREATE TASK\n");
  AliAnalysisTaskCharmBaryonsMC *task = new AliAnalysisTaskCharmBaryonsMC("TaskCharmBaryonsMC");
  mgr->AddTask(task);

  // Create and connect containers for input/output  
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  TString outputfile = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutput1   = mgr->CreateContainer(Form("CharmBaryonHistoMC"),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); // general histos
  mgr->ConnectOutput(task,1,coutput1);
  return task;

}
