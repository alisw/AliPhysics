AliAnalysisTaskQGSep * AddTaskQGSep(Bool_t useMC = kFALSE, Bool_t useAOD = kFALSE)
{
  // create manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("My test train");
  
  // create our task
  AliAnalysisTaskQGSep *task = new AliAnalysisTaskQGSep("AliAnalysisTaskQGSep");

  // uncomment this to use MC
  task->UseMC(useMC);
  task->UseAOD(useAOD);
  
  // create output container
  AliAnalysisDataContainer *output1 = mgr->CreateContainer("cCustomList", TList::Class(), AliAnalysisManager::kOutputContainer,  
							   Form("%s:PWG4_QGSep_UA104",AliAnalysisManager::GetCommonFileName()));

  // add our task to the manager
  mgr->AddTask(task);

  // finaly connect input and output
  mgr->ConnectInput(task, 0,  mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output1);

  return task;

}
