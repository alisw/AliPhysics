void AddTask_GammaJetsPureMC(
                          TString strEffiNeutral = "0.7", 
                          TString strEffiCharged = "0.8", 
                          TString strPDGNonMeas = "-1", 
                          double energyshift = 1.,
                          double jetR = 0.4,                      // jet radius parameter
                          double jetAccEta = 0.8,                 // eta accepatance for jets
                          double jetAccEtaPartPtSpec = 0.4,       // eta acceptance for particles that contribute to the pT spectrum measurement
                          double jetAccEtaPartFF = 1.2,           // eta acceptance for particles that contribute to the FF measurement
                          TString nameTask = "GammaJetsPureMC"
                        ) {

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
  task= new AliAnalysisTaskGammaJetsPureMC(nameTask);

  task->SetEfficiency(strEffiNeutral, strEffiCharged);
  if(!strPDGNonMeas.EqualTo("")) task->SetParticlesMeas(strPDGNonMeas);
  if(energyshift != 1.) task->SetJetEnergyShift(energyshift);
  task->SetJetR(jetR);
  task->SetJetEtaAcc(jetAccEta);
  task->SetEtaAccPartPtSpec(jetAccEtaPartPtSpec);
  task->SetEtaAccPartFF(jetAccEtaPartFF);

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(nameTask, TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), nameTask.Data()));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
