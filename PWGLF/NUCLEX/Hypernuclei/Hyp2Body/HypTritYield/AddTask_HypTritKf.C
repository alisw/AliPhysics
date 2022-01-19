AliAnalysisTask *AddTask_HypTritKf(UInt_t triggerMask = AliVEvent::kINT7 | AliVEvent::kTRD | AliVEvent::kHighMultV0 | AliVEvent::kHighMultSPD, Bool_t pidQa = kFALSE, Bool_t betheSplines = kFALSE, Int_t period = 2017) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_HypTritKf", "No analysis manager found.");
    return 0;
  }
  AliAnalysisTaskHypTritKf *task = new AliAnalysisTaskHypTritKf("mhartungTaskHypTritKf");
  task->SelectCollisionCandidates(triggerMask);
  task->SetTriggerMask(triggerMask);
  task->SetPeriod(period);
  task->SetBetheSplines(betheSplines);
/*
  task->SetUseExternalSplines(kTRUE);
  Double_t he3Params[6];
  Double_t tParams[6];
  // He3
  he3Params[0] = 3.77351;
  he3Params[1] = 12.8008;
  he3Params[2] = 0.00607617;
  he3Params[3] = 2.45706;
  he3Params[4] = 3.21805;
  he3Params[5] = 0.06;
  // Triton
  tParams[0] = 2.51309;
  tParams[1] = 19.2399;
  tParams[2] = 2.89436;
  tParams[3] = 2.09173;
  tParams[4] = -2.2805;
  tParams[5] = 0.06;
  task->SetParamsHe(he3Params);
  task->SetParamsT(tParams);
*/
  //Data Containers
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("histogramsKF", TList::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput2 =
    mgr->CreateContainer("treeKF", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput3 =
    mgr->CreateContainer("treeKF_mc", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  mgr->ConnectOutput(task, 3, coutput3);
  return task;
}
