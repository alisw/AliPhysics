void AddTask_GammaJetsPureMC(TString strEffiNeutral = "0.7", TString strEffiCharged = "0.8", TString strPDGNonMeas = "-1") {

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_GammaJetsPureMC", "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  TString type = mgr->GetInputEventHandler()->GetDataType();
  printf("data type: %s\n", type.Data());
  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  //            find input container
  AliAnalysisTaskGammaJetsPureMC *task=NULL;
  task= new AliAnalysisTaskGammaJetsPureMC("GammaJetsPureMC");

  task->SetEfficiency(strEffiNeutral, strEffiCharged);
  if(!strPDGNonMeas.EqualTo("")) task->SetParticlesNonMeas(strPDGNonMeas);

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaJetsPureMC"), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:GammaJetsPureMC",AliAnalysisManager::GetCommonFileName()));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
