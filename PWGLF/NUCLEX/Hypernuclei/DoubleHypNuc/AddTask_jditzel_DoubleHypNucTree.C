AliAnalysisTask *AddTask_jditzel_DoubleHypNucTree(UInt_t triggerMask = AliVEvent::kAny, Bool_t betheSplines = kFALSE, Int_t period = 2015) {
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
  
  Double_t he3Params[6];
  Double_t tParams[6];

  
  if(period == 2015){
    he3Params[0] = 2.45605;
    he3Params[1] = 19.8067;
    he3Params[2] = -0.77472;
    he3Params[3] = 1.96279;
    he3Params[4] = 0.172695;
    he3Params[5] = 0.05655;
    
    tParams[0] = 2.32603;
    tParams[1] = 19.2492;
    tParams[2] = 30.7943;
    tParams[3] = 2.1697;
    tParams[4] = -8.11114;
    tParams[5] = 0.09311;
  }
  if(period == 2016){
    he3Params[0] = 1.69155;
    he3Params[1] =  27.4992;
    he3Params[2] =  4.00313e-15;
    he3Params[3] =  2.48485;
    he3Params[4] = 8.31768;
    he3Params[5] =  0.05655;
    
    tParams[0] = 1.69155;
    tParams[1] = 27.4992;
    tParams[2] =  4.00313e-15;
    tParams[3] = 2.48485;
    tParams[4] = 8.31768;
    tParams[5] = 0.09311;
  }
  
  if (period == 2017){
    he3Params[0] = 3.20025;
    he3Params[1] = 16.4971;
    he3Params[2] = -0.0116571;
    he3Params[3] = 2.3152;
    he3Params[4] = 3.11135;
    he3Params[5] = 0.06;
    
    tParams[0] = 0.420434;
    tParams[1] = 106.102;
    tParams[2] = -3.15587e-07;
    tParams[3] = 2.32499;
    tParams[4] = 21.3439;
    tParams[5] = 0.06;
  }
  
  if(period == 2018 || period == 2016){
    tParams[0] = 0.427978;
    tParams[1] = 105.46;
    tParams[2] =-7.08642e-07;
    tParams[3] = 2.23332;
    tParams[4] = 18.8231;
    tParams[5] = 0.06;
    
    he3Params[0] = 1.81085;
    he3Params[1] = 29.4656;
    he3Params[2] = 0.0458225;
    he3Params[3] = 2.08689;
    he3Params[4] = 2.28772;
    he3Params[5] = 0.06;
  }
  task->SetParamsHe(he3Params);
  task->SetParamsT(tParams);
  //Data Containers
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histograms", TList::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("tree", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  return task;
}
