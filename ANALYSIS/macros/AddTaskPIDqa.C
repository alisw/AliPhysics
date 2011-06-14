AliAnalysisTask *AddTaskPIDqa(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskPIDqa", "No analysis manager found.");
    return 0;
  }
  //============= Set Task Name ===================
  TString taskName=("AliAnalysisTaskPIDqa");
  //===============================================
  //            Load the task
  gROOT->LoadMacro(Form("%s.cxx+",taskName.Data()));
  
  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskPIDqa *task=new AliAnalysisTaskPIDqa("PIDqaTask");
  mgr->AddTask(task);

  
  //================================================
  //              data containers
  //================================================

  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("PIDqa", TList::Class(),
                         AliAnalysisManager::kOutputContainer,"PIDqa.root");
  
  //           connect containers
  mgr->ConnectInput  (task,  0, mgr->GetCommonInputContainer() );
  mgr->ConnectOutput (task,  1, coutput1);
  
  return task;
}
