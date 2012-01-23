AliAnalysisTaskMuonTreeBuilder *AddTaskTreeBuilder(Bool_t ismc=kFALSE, Int_t run_num=0){
  printf("Inside add task\n");
  // Get the pointer to the existing analysis manager via the static access method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMuonTreeBuilder", "No analysis manager to connect to");
    return NULL;
  }   

  // MC handler if needed
  if(ismc){
  AliMCEventHandler *mcH = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
  if (!mcH) {
    ::Error("AddTaskTreeBuilder", "No MC handler connected");
    return NULL;
  }	
  }

  // The task
  AliAnalysisTaskMuonTreeBuilder *task = new AliAnalysisTaskMuonTreeBuilder("AliAnalysisTaskMuonTreeBuilder");
  if(ismc) task->SetIsMC(kTRUE);

  //outputs
//   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist0",TList::Class(),AliAnalysisManager::kOutputContainer,"final02.root");
//   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("ctree0",TTree::Class(),AliAnalysisManager::kOutputContainer,"final02.root");
  char outname[30];
  sprintf(outname,"TreeRUN%d.root",run_num);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("ctree0",TTree::Class(),AliAnalysisManager::kOutputContainer,outname);

  // Adding the task to the analysis manager
  mgr->AddTask(task);
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
//   mgr->ConnectOutput(task,2,coutput2);

  return task;
}
