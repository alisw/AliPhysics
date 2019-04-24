void AddTask_GammaMCStudies(Double_t maxpT = 100,Int_t nevents = 1) {

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_GammaMCStudies", "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  //            find input container
  AliAnalysisTaskGammaMCStudies *task=NULL;
  task= new AliAnalysisTaskGammaMCStudies("GammaMCStudies");

  task->SetMaxPt(maxpT);
  task->SetNEvents(nevents);
  task->SetNRejectEvents(10000);

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer("GammaStudies", TList::Class(), AliAnalysisManager::kOutputContainer,AliAnalysisManager::GetCommonFileName());

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
