//________________________________________________________________________
void demoBatch() {
  const char *collectionfile = "wn.xml";

  //____________________________________________//
  //Usage of event tags
  AliTagAnalysis *analysis = new AliTagAnalysis();
  TChain *chain = 0x0;
  chain = analysis->GetInputChain("pp",collectionfile);

  //____________________________________________//
  //Usage of a collection file - no tags
  /*TChain *chain = new TChain("esdTree");
  TAlienCollection *collection = TAlienCollection::Open(collectionfile);
  TGridResult *result = collection->GetGridResult("",0,0);  
  for(Int_t i = 0; i < resukt->GetEntries(); i++)  chain->Add(result->GetKey(i,"turl"));*/
  
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
  
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("grid",chain);
  }
}                         
                      
