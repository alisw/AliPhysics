AliAnalysisTask* AddTaskCL(Bool_t mcmode=0, Double_t ptmin=3, Double_t m2cut=2, Bool_t doSkip=0, Bool_t doClus=1)
{  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType();

  TString name("CLTask");
  AliCLTask *task = new AliCLTask(name);
  task->SetMCMode(mcmode);
  task->SetPtCut(ptmin);
  task->SetM2Cut(m2cut);
  task->SetDoSkip(doSkip);
  task->SetDoClust(doClus);

  mgr->AddTask(task);
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer()  ;
  TString contName(name);
  contName += "_histos";
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(contName.Data(), 
							   TTree::Class(),
							   AliAnalysisManager::kOutputContainer,
							   Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (task, 0,  cinput );
  mgr->ConnectOutput (task, 1, coutput );
  return task;
}
