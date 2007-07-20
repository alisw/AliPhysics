//________________________________________________________________________
void demoInteractive() {
  //____________________________________________//
  AliTagAnalysis *TagAna = new AliTagAnalysis();
 
  AliRunTagCuts *runCuts = new AliRunTagCuts();
  AliLHCTagCuts *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts *evCuts = new AliEventTagCuts();
  evCuts->SetMultiplicityRange(11,12);

  //grid tags
  TAlienCollection* coll = TAlienCollection::Open("../tag.xml");
  TGridResult* TagResult = coll->GetGridResult("",0,0);
  TagAna->ChainGridTags(TagResult);
  TChain* chain = 0x0;
  chain = TagAna->QueryTags(runCuts,lhcCuts,detCuts,evCuts);
 
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  //____________________________________________//
  // 1st Pt task
  AliAnalysisTaskPt *task1 = new AliAnalysisTaskPt("TaskPt");
  mgr->AddTask(task1);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist1", TH1::Class(),AliAnalysisManager::kOutputContainer,"Pt.ESD.root");
  
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  cinput1->SetData(chain);
  
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain);
  }
}                         
                      
