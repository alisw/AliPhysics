void run(const Char_t *list=0x0, Int_t nmax=-1) {
  TStopwatch timer;
  timer.Start();

  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libTRDqaRec.so");
	

  //____________________________________________//

  gROOT->LoadMacro(Form("%s/TRD/qaRec/CreateESDChain.C", gSystem->ExpandPathName("$ALICE_ROOT")));
  TChain *chain = CreateESDChain(list, nmax);
  //chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("*FMD*",0);
  chain->SetBranchStatus("*Calo*",0);
  chain->SetBranchStatus("Tracks", 1);
  chain->SetBranchStatus("ESDfriend*",1);
  chain->Lookup();
  chain->GetListOfFiles()->Print();
  printf("\n ----> CHAIN HAS %d ENTRIES <----\n\n", (Int_t)chain->GetEntries());
  
  AliLog::SetGlobalLogLevel(AliLog::kError);

  //____________________________________________
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TRD Track Info Manager");
  //mgr->SetSpecialOutputLocation(source); // To Be Changed
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);  
  AliMCEventHandler *mc = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mc);
  //mgr->SetDebugLevel(10);

  //____________________________________________
  // TRD track summary generator
  AliTRDtrackInfoGen *task1 = new AliTRDtrackInfoGen();
  task1->SetDebugLevel(1);
  mgr->AddTask(task1);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("data", TChain::Class(), AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("TrackInfoList", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  mgr->ConnectInput( task1, 0, cinput1);
  mgr->ConnectOutput(task1, 0, coutput1);

  //____________________________________________
  // TRD barrel tracking efficiency
  AliTRDtrackingEfficiency *task2 = new AliTRDtrackingEfficiency();
  task2->SetDebugLevel(1);
  mgr->AddTask(task2);
  //Create containers for input/output
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("TrackingEfficiency", TList::Class(), AliAnalysisManager::kOutputContainer, "TRD.TrackingEfficiency.root");
  mgr->ConnectInput( task2, 0, coutput1);
  mgr->ConnectOutput(task2, 0, coutput2);

  //____________________________________________
  // TRD combined tracking efficiency
  AliTRDtrackingEfficiencyCombined *task3 = new AliTRDtrackingEfficiencyCombined();
  task3->SetDebugLevel(0);
  mgr->AddTask(task3);
  // Create containers for input/output
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("TrackingEfficiencyCombined", TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.TrackingEfficiencyCombined.root");
  mgr->ConnectInput( task3, 0, coutput1);
  mgr->ConnectOutput(task3, 0, coutput3);



  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}
