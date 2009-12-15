

void AddSDDPoints(Int_t run){
  gROOT->LoadMacro("AliAnalysisTaskSDDRP.cxx+g");
  AliAnalysisTaskSDDRP *task= new AliAnalysisTaskSDDRP();
  task->SetGeometryFile("geometry.root");
  task->SetRunNumber(run);
  mgr->AddTask(task);
  mgr->SetDebugLevel(2);


  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("coutputRP",TList::Class(),AliAnalysisManager::kOutputContainer,"/home/prino/alice/test/SDDPoints.root");



  mgr->ConnectInput(task,0,cinput1);
  mgr->ConnectOutput(task,0,coutput1);


  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", esdChain);
  }   
  delete mgr;
}   



