AliAnalysisTask *AddTask_jditzel_DoubleHypNucTree(UInt_t triggerMask = AliVEvent::kAny, Bool_t betheSplines = kFALSE, Int_t period = 2015, Int_t methodnum = 1) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_jditzel_DoubleHypNucTree", "No analysis manager found.");
    return 0;
  }
  
  AliAnalysisTaskDoubleHypNucTree *task = new AliAnalysisTaskDoubleHypNucTree("jditzelTaskDoubleHypNucTree");
  
  task->SelectCollisionCandidates(triggerMask);
  task->SetTriggerMask(triggerMask);
  task->SetPeriod(period);
  task->SetBetheSplines(betheSplines);
  task->SetMethod(methodnum);
  
  Double_t he3Params[6];
  Double_t tParams[6];

  
  if(period == 2015){
    he3Params[0] = 2.45605;
    he3Params[1] = 19.8067;
    he3Params[2] = -0.77472;
    he3Params[3] = 1.96279;
    he3Params[4] = 0.172695;
    he3Params[5] = 0.06;
    
    tParams[0] = 2.32603;
    tParams[1] = 19.2492;
    tParams[2] = 30.7943;
    tParams[3] = 2.1697;
    tParams[4] = -8.11114;
    tParams[5] = 0.06;
  }
  
  if(period == 2018 || period == 2016){
    tParams[0] = 0.669634;
    tParams[1] = 53.1497;
    tParams[2] =-1.32853e-08;
    tParams[3] = 2.5775;
    tParams[4] = 17.7607;
    tParams[5] = 0.06;
    
    he3Params[0] = 1.50582;
    he3Params[1] = 33.7232;
    he3Params[2] = -0.0923749;
    he3Params[3] = 2.00901;
    he3Params[4] = 2.28772;
    he3Params[5] = 0.06;
  }
  task->SetParamsHe(he3Params);
  task->SetParamsT(tParams);
  //Data Containers
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histograms", TList::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("tree4LH", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("tree4Li", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("tree5LHe", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("tree5LHe2", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput6 = mgr->CreateContainer("tree4LHe", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput7 = mgr->CreateContainer("tree4LH3B", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput8 = mgr->CreateContainer("tree4LLH", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  mgr->ConnectOutput(task, 3, coutput3);
  mgr->ConnectOutput(task, 4, coutput4);
  mgr->ConnectOutput(task, 5, coutput5);
  mgr->ConnectOutput(task, 6, coutput6);
  mgr->ConnectOutput(task, 7, coutput7);
  mgr->ConnectOutput(task, 8, coutput8);
  return task;
}
