AliAnalysisTask *AddTaskPIDqa(const char *useroutputfile=""){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskPIDqa", "No analysis manager found.");
    return 0;
  }
  
  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskPIDqa *task=new AliAnalysisTaskPIDqa("PIDqaTask");
  mgr->AddTask(task);
  
  //================================================
  //              data containers
  //================================================

  TString outputfile=useroutputfile;
  if (outputfile.IsNull()) outputfile = Form("%s:PIDqa", AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("PIDqa", TList::Class(),
                         AliAnalysisManager::kOutputContainer,outputfile);
  
  //           connect containers
  mgr->ConnectInput  (task,  0, mgr->GetCommonInputContainer() );
  mgr->ConnectOutput (task,  1, coutput1);
  
  return task;
}
