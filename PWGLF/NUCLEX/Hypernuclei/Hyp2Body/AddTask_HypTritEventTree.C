AliAnalysisTask *AddTask_HypTritEventTree(UInt_t triggerMask = AliVEvent::kINT7 | AliVEvent::kTRD | AliVEvent::kHighMultV0 | AliVEvent::kHighMultSPD, Bool_t pidQa = kFALSE, Bool_t betheSplines = kFALSE, Int_t period = 2017) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_HypTritEventTree", "No analysis manager found.");
    return 0;
  }
  AliAnalysisTaskHypTritEventTree *task = new AliAnalysisTaskHypTritEventTree("mhartungTaskHypTritEventTree");
  task->SelectCollisionCandidates(triggerMask);
  task->SetTriggerMask(triggerMask);
  task->SetPidQa(pidQa);
  task->SetPeriod(period);
  task->SetBetheSplines(betheSplines);
  if(period == 2016){
    Double_t he3Params[6] = {0.715489, 59.5463, 4.44487e-12, 2.69874, 24.063, 0.04725};
    Double_t tParams[6] = {0.223948, 180.564, -3.03884e-10, 2.30095, 34.2269, 0.06517};
  }
  task->SetParamsHe(he3Params);
  task->SetParamsT(tParams);
  //Data Containers
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("histograms", TList::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput2 =
    mgr->CreateContainer("tree3LH", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput3 =
    mgr->CreateContainer("tree3LHGen", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  mgr->ConnectOutput(task, 3, coutput3);
  return task;
}
