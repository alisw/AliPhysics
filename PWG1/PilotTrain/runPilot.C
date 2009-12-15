void runPilot(Int_t run) {
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
  
 

  

  Bool_t doQAsym        = 1;
  Bool_t doVZERO        = 1;
  Bool_t doVertex       = 1;
  Bool_t doSPD          = 1;  
  Bool_t doFMD          = 1;
  Bool_t doTPC          = 1;
  Bool_t doEventStat    = 1;
  Bool_t doSDD          = 1;
   //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  mgr->SetDebugLevel(2);
  
  AliVEventHandler* esdH = new AliESDInputHandlerRP;
  mgr->SetInputEventHandler(esdH);  

  // Chain 
  TChain* chain = new TChain("esdTree");
  chain->AddFile("~/104321/AliESDs.root");
  
  //
  // Wagons
  //
  //
  // Collision Selector (static)
  AliPhysicsSelection* colsel =  new AliPhysicsSelection();
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
      task2->SelectCollisionCandidates();
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
      gROOT->LoadMacro("AliAnalysisTaskSDDRP.cxx++g");
      gROOT->LoadMacro("AddSDDPoints.C");
      AliAnalysisTaskSE* task5 = AddSDDPoints(run);
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
      AddTaskEventStats();
  }

  
  
  // Init
  if (!mgr->InitAnalysis()) 
      mgr->PrintStatus();
  
  // Run on dataset
  mgr->StartAnalysis("local", chain, 500);
  timer.Stop();
  timer.Print();
}

