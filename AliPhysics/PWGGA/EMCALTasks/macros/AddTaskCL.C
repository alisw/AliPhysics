AliAnalysisTask* AddTaskCL(Bool_t mcmode=0, Double_t ptmin=3, Double_t m2cut=2, Bool_t doSkip=0, Bool_t doClus=1, Bool_t doList=0)
{  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType();

  TString name("CLTask");
  AliCLTask *task = new AliCLTask(name,doList);
  task->SetMCMode(mcmode);
  task->SetPtCut(ptmin);
  task->SetM2Cut(m2cut);
  task->SetDoSkip(doSkip);
  task->SetDoClust(doClus);
  mgr->AddTask(task);
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer()  ;
  TString contName1(name);
  contName1 += "_tree";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), 
							    TTree::Class(),
							    AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (task, 0, cinput);
  mgr->ConnectOutput (task, 1, coutput1);
  if (doList) {
    TString contName2(name);
    contName2 += "_hists";
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(), 
							      TList::Class(),
							      AliAnalysisManager::kOutputContainer,
							      Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput (task, 2, coutput2);
  }
  return task;
}
