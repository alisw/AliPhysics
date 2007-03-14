//________________________________________________________________________
void demoLocal() {
  TChain* chain = new TChain("esdTree");
  chain->Add("data/001/AliESDs.root");
  chain->Add("data/002/AliESDs.root");

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  //____________________________________________//
  // 1st Pt task
  AliAnalysisTaskRLPt *task1 = new AliAnalysisTaskRLPt("TaskRLPt");
  mgr->AddTask(task1);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist1", TH1::Class(),AliAnalysisManager::kOutputContainer,"Pt.ESD.root");
  
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  cinput1->SetData(chain);
  
  if(mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain);
  }
}                         
                      
