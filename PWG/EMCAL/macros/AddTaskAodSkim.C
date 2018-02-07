// $Id: AddTaskSkim.C 4586 2013-01-16 15:32:16Z loizides $

AliAodSkimTask *AddTaskAodSkim(Double_t mine=5) 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskAodSkim", "No analysis manager found.");
    return 0;
  }

  AliAODHandler *output = new AliAODHandler;
  output->SetOutputFileName("aodskim.root");
  output->SetFillAOD(1);
  output->SetFillAODforRun(1);
  output->SetFillExtension(0);
  mgr->SetOutputEventHandler(output);

  const char *name = "skimtask";
  AliAodSkimTask *task = new AliAodSkimTask(name);
  task->SetClusMinE(mine);
  
  mgr->AddTask(task);
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer()  ;
  TString contName(name);
  contName += "_histos";
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(contName.Data(), 
							   TList::Class(),
							   AliAnalysisManager::kOutputContainer,
							   Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (task, 0,  cinput );
  mgr->ConnectOutput (task, 1, coutput );

  return task;
}
