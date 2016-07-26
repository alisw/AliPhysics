AliAnalysisTask *AddTaskTRDPIDmatching(Int_t trigger=0, Int_t system=0){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_pachmay_PIDTRDTree", "No analysis manager found.");
    return 0;
  }


  //========= Add task to the ANALYSIS manager =====
  AliTRDPIDmatching *task = new AliTRDPIDmatching("TRDPIDmatching");
  if(trigger==0) task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kINT1 | AliVEvent::kINT5 | AliVEvent::kINT7 | AliVEvent::kINT8);
  if(trigger==1) task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kINT7 | AliVEvent::kINT8);
  if(trigger==2) task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  if(trigger==3) task->SelectCollisionCandidates(AliVEvent::kTRD);
  if(system==1) task->SetppAnalysis();
  if(system==2) task->SetpPbAnalysis();
  if(system==3) task->SetPbPbAnalysis();
  mgr->AddTask(task);


  //================================================
  //              data containers
  //================================================
  //            find input container
  //below the trunk version
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  AliAnalysisDataContainer *coutput1 =
      mgr->CreateContainer("TRDPIDmatching",
                           TList::Class(),
                           AliAnalysisManager::kOutputContainer,
                           "TRDPIDmatching.root");

  //connect containers
  mgr->ConnectInput(task,  0, cinput  );
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
}
