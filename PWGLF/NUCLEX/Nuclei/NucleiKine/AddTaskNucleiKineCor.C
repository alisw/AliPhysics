AliAnalysisTaskNucleiKineCor* AddTaskNucleiKineCor(Double_t trigpt=5, TString contname="cor")
{
  AliAnalysisTaskNucleiKineCor *task = new AliAnalysisTaskNucleiKineCor();
  task->SetPt(trigpt);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    cout<<"AliAnalysisTaskNucleiKine: No analysis manager to connect to."<<endl;
    return NULL;
  }
  task->SetName(Form("task_%s_%.1f",contname.Data(),trigpt));
  mgr->AddTask(task);

  // Create containers for input/output
  TString outputFileName = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("output_%s_%.1f",contname.Data(),trigpt), 
							   TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);

  //connect containers
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  return task;
}
