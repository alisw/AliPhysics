AliAnalysisTask *AddTask_doenigus_HdibaryonLPpi(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_doenigus_HdibaryonLPpi", "No analysis manager found.");
    return 0;
  }

  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskSE *taskHdibaryonLPpi = new AliAnalysisTaskHdibaryonLPpi("doenigus_HdibaryonLPpi");

  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  mgr->AddTask(taskHdibaryonLPpi);

  if (!hasMC){
    taskHdibaryonLPpi->SelectCollisionCandidates(AliVEvent::kAny);
  }

  //================================================
  //              data containers
  //================================================
  //            find input container

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  AliAnalysisDataContainer *coutput_doenigusHdibaryonLPpi = mgr->CreateContainer("doenigus_HdibaryonLPpi", TList::Class(), AliAnalysisManager::kOutputContainer,"doenigus_HdibaryonLPpi.root");

  //           connect containers
  mgr->ConnectInput  (taskHdibaryonLPpi,  0, cinput );
  mgr->ConnectOutput (taskHdibaryonLPpi,  1, coutput_doenigusHdibaryonLPpi);

  return taskHdibaryonLPpi;
}
