void AddTask_EtaPrimeMCStudies(TString specialSetting = "EtaPrimeBiased") {

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_EtaPrimeMCStudies", "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  //            find input container
  AliAnalysisTaskEtaPrimeMCStudies *task = NULL;
  task = new AliAnalysisTaskEtaPrimeMCStudies("EtaPrimeMCStudies");
  if( !task ){
    Error("AddTask_EtaPrimeMCStudies","Failed to create the task.");
    return;
  }

  //connect containers
  AliAnalysisDataContainer *coutput = mgr->CreateContainer( "EtaPrimeMCStudies", TList::Class(), AliAnalysisManager::kOutputContainer, Form( "EtaPrimeMCStudies_%s.root", specialSetting.Data()));
  AliAnalysisDataContainer *couttree = mgr->CreateContainer( "EtaPrimeMCStudiesTree", TTree::Class(), AliAnalysisManager::kOutputContainer, Form( "EtaPrimeMCStudiesTree_%s.root", specialSetting.Data() ));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  mgr->ConnectOutput(task,2,couttree);

  return;

}
