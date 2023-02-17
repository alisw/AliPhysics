AliAnalysisTaskNFMs* AddTaskNFMs(TString name, UInt_t fbit, int mincnt, int maxcnt, float pt0, float pt1, float pt2, float pt3, float pt4, float pt5, float pt6, float pt7, bool fmmax, bool fpileup, bool fncl, bool fncr, float vncl, float vncr, bool run1)
{
  // get the manager via the static access member. since it's static, you don't need
  // an instance of the class to call the function
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }   
  // get the input event handler, again via a static method. 
  // this handler is part of the managing system and feeds events
  // to your task
  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }

  // by default, a file is open for writing. here, we get the filename
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":Intermittency";      // create a subfolder in the file

  // now we create an instance for task
  AliAnalysisTaskNFMs *task = new AliAnalysisTaskNFMs(name.Data());
  if(!task) return 0x0;
  task->SetFilterBit(fbit);
  task->SetCentLim(mincnt, maxcnt);                     //min and max cent value
  task->SetpTbins(pt0, pt1, pt2, pt3, pt4, pt5, pt6, pt7);  //pt bin vals in order
  task->SetMmax(fmmax);                                   //kTRUE for m upto 82
  task->SetPileUpRead(fpileup);                           //pileup flag
  task->SetSysflag(fncl, fncr, vncl, vncr);               //ncluster, ncrossdrows flag, val
  task->SetDataset(run1);                                 //kTRUE for Run1
  printf("Container name is %s\n",name.Data());

  // add your task to the manager
  mgr->AddTask(task);
  // your task needs input: here we connect the manager to your task
  AliAnalysisDataContainer  *cinput = mgr->GetCommonInputContainer();
  // same for the output 
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("outputlist", TList::Class(),AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("ntplist1", TList::Class(),AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("ntplist2", TList::Class(),AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("ntplist3", TList::Class(),AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("ntplist4", TList::Class(),AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput6 = mgr->CreateContainer("QAlist", TList::Class(),AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());

  //connect containers
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  mgr->ConnectOutput(task, 3, coutput3);
  mgr->ConnectOutput(task, 4, coutput4);
  mgr->ConnectOutput(task, 5, coutput5);
  mgr->ConnectOutput(task, 6, coutput6);
  // in the end, this macro returns a pointer to your task. this will be convenient later on
  // when you will run your analysis in an analysis train on grid
  return task;
}
