//________________________________________________________________________
void demoBatch(const char* collectionfile) {
  //____________________________________________//
  // Open a collection
  TAlienCollection* coll = TAlienCollection::Open(collectionfile);
  if (!coll) Error("demoBatch", Form("Cannot create an AliEn collection from %s", collectionfile));

  //____________________________________________//
  // Create an esd chain
  const char *chainname="esdTree";
  TChain* chain1 = new TChain(chainname);
  
  //____________________________________________//
  // Convert the collection and fill tha chain
  coll->Reset();
  while ( coll->Next() ) {
    Info("demoBatch", Form("Adding %s", coll->GetTURL("")));
    chain1->Add(coll->GetTURL(""));
  }

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
                      
