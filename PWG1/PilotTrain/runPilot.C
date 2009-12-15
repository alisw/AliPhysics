void runPilot() {
  TStopwatch timer;
  timer.Start();
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libPWG0base.so");
  gSystem->Load("libPWG0dep.so");
  gSystem->Load("libPWG0selectors.so");
  gSystem->Load("libPWG1.so");
  gSystem->Load("libPWG2.so");
  gSystem->Load("libPWG2forward.so");
  
 

  gROOT->LoadMacro("AliESDInputHandlerRPITS.cxx++g");  

  Bool_t doQAsym        = 1;   // output ok
  Bool_t doVZERO        = 1;   // output ok but there is a 2nd file
  Bool_t doVertex       = 1;   // output ok
  Bool_t doSPD          = 1;   // output ok, needs RP   
  Bool_t doFMD          = 1;   // output ok
  Bool_t doTPC          = 1;   // output ok
  Bool_t doEventStat    = 1;   // output ok
  Bool_t doSDD          = 1;   // outout ok needs RP
   //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  mgr->SetDebugLevel(2);
  

  AliInputEventHandler* esdH = new AliESDInputHandlerRPITS();
  esdH->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdH);  

  // Chain 
  TGrid::Connect("alien://");
  
  TChain* chain = new TChain("esdTree");
  chain->AddFile("alien:///alice/data/2009/LHC09d/000104321/ESDs/pass1/09000104321018.10/AliESDs.root");
  chain->AddFile("alien:///alice/data/2009/LHC09d/000104321/ESDs/pass1/09000104321018.20/AliESDs.root");
  //
  // Wagons
  //
  //
  // Collision Selector (static)
  AliPhysicsSelection* colsel =  new AliPhysicsSelection();
  colsel->AddBackgroundIdentification(new AliBackgroundSelection());

  AliAnalysisTaskSE::SetCollisionSelector(colsel);
  

  // TPC QA (E. Sicking)
  //
  if (doQAsym) {
      gROOT->LoadMacro("AddTaskQAsym.C");
      AliAnalysisTaskSE* task1 = AddTaskQAsym();
      task1->SelectCollisionCandidates();
  }
  
  //
  // VZERO QA  (C. Cheshkov)
  //
  if (doVZERO) {
      gROOT->LoadMacro("AddTaskVZEROQA.C");
      AliAnalysisTaskSE* task2 =  AddTaskVZEROQA(0);
//      task2->SelectCollisionCandidates();
  }
  
  //
  // Vertexing (A. Dainese)
  // 
  if (doVertex) {
      gROOT->LoadMacro("$(ALICE_ROOT)/PWG1/macros/AddTaskVertexESD.C");
      AliAnalysisTaskVertexESD* task3 =  AddTaskVertexESD();
      task3->SelectCollisionCandidates();
  }
  
  //
  // SPD (A. Mastroserio)
  //
  if (doSPD) {
      gROOT->LoadMacro("AliAnalysisTaskSPD.cxx++g");
      gROOT->LoadMacro("AddTaskSPDQA.C");
      AliAnalysisTaskSE* task4 = AddTaskSPDQA();
      task4->SelectCollisionCandidates();
  }

  //
  // SDD (F. Prino)
  //
  if (doSDD) {
      gROOT->LoadMacro("AddSDDPoints.C");
      AliAnalysisTaskSE* task5 = AddSDDPoints();
      task5->SelectCollisionCandidates();
  }
  
  //
  // FMD (Hans Hjersing Dalsgaard)
  //
  if (doFMD) {
      gROOT->LoadMacro("AddTaskFMD.C");
      AliAnalysisTaskSE* task6 = AddTaskFMD();
      task6->SelectCollisionCandidates();
  }
  
  //
  // TPC (Jacek Otwinowski)
  //
  if (doTPC) {
      gROOT->LoadMacro("$(ALICE_ROOT)/PWG1/TPC/macros/AddTaskPerformanceTPCQA.C");
      AliPerformanceTask *tpcQA = AddTaskPerformanceTPCQA(kFALSE, kTRUE);
  }

  //
  // Event Statistics (Jan Fiete)
  //

  if (doEventStat) {
      gROOT->LoadMacro("AddTaskEventStats.C");
      evtStats = AddTaskEventStats();
      evtStats->SetPhysicsSelection(colsel);
      AliLog::SetClassDebugLevel("AliPhysicsSelection", AliLog::kDebug);
  }

  
  
  // Init
  if (!mgr->InitAnalysis()) 
      mgr->PrintStatus();
  
  // Run on dataset
  mgr->StartAnalysis("local", chain);
  timer.Stop();
  timer.Print();
}

