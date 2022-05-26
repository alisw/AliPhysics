AliAnalysisTaskMSDibaryons *AddTaskMSDibaryons(TString taskname="MSDibaryons")
{

  // Connection to the analysis manager
  //================================================================
  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    printf("ERROR: No analysis manager to connect to.\n");
    return NULL;
  }
  
  // Task creation
  //================================================================
  AliAnalysisTaskMSDibaryons *task=new AliAnalysisTaskMSDibaryons(taskname.Data());
  UInt_t trigSel=AliVEvent::kINT7;
  task->SetTrigger(trigSel);
  mgr->AddTask(task);

  // Connect input/output
  //================================================================
  AliAnalysisDataContainer *cinput =mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput=mgr->CreateContainer("cList",TList::Class(),
							 AliAnalysisManager::kOutputContainer,
							 AliAnalysisManager::GetCommonFileName());
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return task;

}
