//________________________________________________________________________
void demoCAF(TChain *chain, const char* mode) {
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  //____________________________________________//
  // 1st Pt task
  AliAnalysisTask *task1 = new AliAnalysisTaskPt("TaskPt");
  mgr->AddTask(task1);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist1", TH1::Class(),AliAnalysisManager::kOutputContainer, "Pt.ESD.1.root");
  //____________________________________________//
                                                                                                                                              
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  if (mgr->InitAnalysis()) mgr->PrintStatus();
  else return;

  mgr->StartAnalysis(mode,chain);
}                         
                      
