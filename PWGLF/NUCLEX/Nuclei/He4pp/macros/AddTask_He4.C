AliAnalysisTask *AddTask_He4(){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_He4", "No analysis manager found.");
    return 0;
  }


  //========= Add task to the ANALYSIS manager =====
  AliAnalysisHe4 *task = new AliAnalysisHe4("TaskHe4");
  task->Initialize();
  mgr->AddTask(task);
  
  
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  //dumm output container
  AliAnalysisDataContainer *coutput = 
      mgr->CreateContainer("List_He4", 
                           TList::Class(),
                           AliAnalysisManager::kOutputContainer,
                           AliAnalysisManager::GetCommonFileName());


  //connect containers
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  1, coutput);


  return task;
}
