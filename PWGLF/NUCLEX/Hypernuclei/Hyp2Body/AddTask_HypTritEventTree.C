AliAnalysisTask *AddTask_HypTritEventTree(UInt_t triggerMask = AliVEvent::kINT7, Bool_t pidQa = kTRUE, Bool_t betheSplines = kTRUE, Int_t period = 2015) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_lkreis_HypTritEventTree", "No analysis manager found.");
    return 0;
  }
  AliAnalysisTaskHypTritEventTree *task = new AliAnalysisTaskHypTritEventTree("lkreisTaskHypTritEventTree");
  task->SelectCollisionCandidates(triggerMask);
  task->SetPidQa(pidQa);
  task->SetPeriod(period);
  task->SetBetheSplines(betheSplines);
  // MC LHC12d3
//  Double_t he3Params[6] = {1.79043,18.0066,0.00212257,2.24677,4.25945,0.06892};
//  Double_t tParams[6] = {2.32603,19.2492,30.7943,2.16971,-8.11114, 0.09311};
  // DATA LHC15o
  Double_t he3Params[6] = {2.45605,19.8067,-0.774729,1.96279,0.172695, 0.05655};
  Double_t tParams[6] = {2.32603,19.2492,30.7943,2.16971,-8.11114, 0.09311};
  task->SetParamsHe(he3Params);
  task->SetParamsT(tParams);
  //Data Containers
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("histograms", TList::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput2 =
    mgr->CreateContainer("tree", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput3 =
    mgr->CreateContainer("tree_mc", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  mgr->ConnectOutput(task, 3, coutput3);
  return task;
}
