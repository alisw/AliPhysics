//________________________________________________________________________
void demoLocal() {
  //____________________________________________//
  AliTagAnalysis *TagAna = new AliTagAnalysis(); 
  TagAna->ChainLocalTags("../Tags/");

  AliEventTagCuts *EvCuts1 = new AliEventTagCuts();
  EvCuts1->SetMultiplicityRange(11,12);  
  TChain* chain1 = 0x0;
  chain1 = TagAna->QueryTags(EvCuts1);

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager();
  //____________________________________________//
  // 1st Pt task
  AliAnalysisTask *task1 = new AliAnalysisTaskPt("TaskPt");
  mgr->AddTask(task1);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist1", TH1::Class(),AliAnalysisManager::kOutputContainer);
  
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  cinput1->SetData(chain1);
  
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    chain1->Process(mgr);
  }
}                         
                      
