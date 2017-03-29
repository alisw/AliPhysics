AliAnalysisTask *AddTaskDeuteronpA(){


  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_janielsk_DeuteronpA", "No analysis manager found.");
    return 0;
  }
  //============= Set Task Name ===================
  //TString taskName=("AliAnalysisDeuteronpA.cxx+g");
  //===============================================
  //            Load the task
  //gROOT->LoadMacro(taskName.Data());


  
  //========= Add task to the ANALYSIS manager =====

  //normal tracks
  AliAnalysisDeuteronpA *task = new AliAnalysisDeuteronpA("janielskTaskDeuteronpA");
  task->SelectCollisionCandidates(AliVEvent::kINT7);

  //initialize task
  task->Initialize();


  //add task to manager
  mgr->AddTask(task);


  

  //================================================
  //              data containers
  //================================================
  //            find input container
  //below the trunk version
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  AliAnalysisDataContainer *coutput1;


  coutput1 =  mgr->CreateContainer(Form("akalweit_DeuteronpA"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:akalweit_DeuteronpA", AliAnalysisManager::GetCommonFileName())); 



  //connect containers

  //
  mgr->ConnectInput  (task,  0, cinput );
  //mgr->ConnectOutput (task,  0, coutput0);
  mgr->ConnectOutput (task,  1, coutput1);

  return task;
}
