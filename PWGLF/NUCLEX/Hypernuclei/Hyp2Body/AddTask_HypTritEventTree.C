AliAnalysisTask *AddTask_HypTritEventTree(UInt_t triggerMask = AliVEvent::kAny, Bool_t pidQa = kFALSE, Bool_t betheSplines = kFALSE, Int_t period = 2017) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_lkreis_HypTritEventTree", "No analysis manager found.");
    return 0;
  }
  AliAnalysisTaskHypTritEventTree *task = new AliAnalysisTaskHypTritEventTree("lkreisTaskHypTritEventTree");
  task->SelectCollisionCandidates(triggerMask);
  task->SetTriggerMask(triggerMask);
  task->SetPidQa(pidQa);
  task->SetPeriod(period);
  task->SetBetheSplines(betheSplines);

/* 
  // MC LHC12d3
//  Double_t he3Params[6] = {1.79043,18.0066,0.00212257,2.24677,4.25945,0.06892};
//  Double_t tParams[6] = {2.32603,19.2492,30.7943,2.16971,-8.11114, 0.09311};
// DATA LHC15o
    if(period == 2015){
    Double_t he3Params[6] = {2.45605,19.8067,-0.774729,1.96279,0.172695, 0.05655};
    Double_t tParams[6] = {2.32603,19.2492,30.7943,2.16971,-8.11114, 0.09311};
  }
    if(period == 2016){
    Double_t he3Params[6] = {1.69155, 27.4992, 4.00313e-15, 2.48485, 8.31768, 0.05655};
    Double_t tParams[6] = {1.69155, 27.4992, 4.00313e-15, 2.48485, 8.31768, 0.09311};
*/ 
    Double_t he3Params[6];
    Double_t tParams[6];

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
