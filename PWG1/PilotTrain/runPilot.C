void runPilot() {
  TStopwatch timer;
  timer.Start();
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT -I$ALICE_ROOT/TRD");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG0base.so");
  gSystem->Load("libPWG0dep.so");
  gSystem->Load("libPWG0selectors.so");
  gSystem->Load("libPWG1.so");
  gSystem->Load("libPWG2.so");
  gSystem->Load("libPWG3muon.so");
  gSystem->Load("libPWG3muondep.so");
  gSystem->Load("libPWG2forward.so");
  gSystem->Load("libPWG4PartCorrBase.so");
  gSystem->Load("libPWG4PartCorrDep.so");
 


  Bool_t doQAsym        = 0;   // (data)
  Bool_t doVZERO        = 0;   // (data)
  Bool_t doVertex       = 0;   // (data)
  Bool_t doSPD          = 0;   // (data/RP)   
  Bool_t doFMD          = 0;   // (data)
  Bool_t doTPC          = 0;   // AliPerformanceTask (data/MCtruth)
  Bool_t doEventStat    = 0;   // (data)
  Bool_t doSDD          = 0;   // (data/RP)
  Bool_t doTRD          = 0;   // TRD 
  Bool_t doITS          = 0;   // ITS
  Bool_t doCALO         = 0;   // Calorimeter
  Bool_t doMUONTrig     = 0;   // MUON trigger
  Bool_t doMUONEff      = 0;   // MUON efficiency  NEEDS geometry
  Bool_t doV0           = 1;   // V0 recosntruction performance NEEDS MCtruth
   
  //
  //
  // Tasks that need MC information
  //

  
  
   //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  mgr->SetDebugLevel(10);
  

  AliInputEventHandler* esdH = new AliESDInputHandlerRP();
  esdH->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdH);  

  // Chain 
  // TGrid::Connect("alien://");
  
  TChain* chain = new TChain("esdTree");
  //chain->AddFile("alien:///alice/data/2009/LHC09d/000104321/ESDs/pass1/09000104321018.10/AliESDs.root");
  //chain->AddFile("alien:///alice/data/2009/LHC09d/000104321/ESDs/pass1/09000104321018.20/AliESDs.root");
  chain->AddFile("./data/AliESDs.root");
  

  //
  // Wagons running on data
  //
  //
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
  // ITS
  // 
  if (doITS) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskPerformanceITS.C");
      AliAnalysisTaskITSTrackingCheck *itsQA = AddTaskPerformanceITS(kFALSE);
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

  // TPC QA (E. Sicking)
  //
  if (doQAsym) {
      gROOT->LoadMacro("AddTaskQAsym.C");
      AliAnalysisTaskSE* task1 = AddTaskQAsym();
      task1->SelectCollisionCandidates();
  }
  //
  // TPC (Jacek Otwinowski)
  //
  if (doTPC) {
      // 
      // Optionally MC information can be used by setting the 1st argument to true
      gROOT->LoadMacro("$(ALICE_ROOT)/PWG1/TPC/macros/AddTaskPerformanceTPCQA.C");
      AliPerformanceTask *tpcQA = AddTaskPerformanceTPCQA(kFALSE, kTRUE);
  }
  
  //
  // TRD (Alex Bercuci, M. Fasel) 
  //
  if(doTRD) {
      AliAnalysisDataContainer *ci[] = {0x0, 0x0, 0x0};
      //
      // Check ESD
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/TRD/macros/AddTRDcheckESD.C++");
      AddTRDcheckESD(mgr);
      //
      // Info top task
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/TRD/macros/AddTRDinfoGen.C++");
      AddTRDinfoGen(mgr, "ALL", 0x0, ci);
      //
      // check DET
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/TRD/macros/AddTRDcheckDET.C++");
      AddTRDcheckDET(mgr, "ALL", ci);
      //
      // check PID (ref maker ???)
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/TRD/macros/AddTRDcheckPID.C++");
      AddTRDcheckPID(mgr, "ALL", ci);
      //
      // Efficiency
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/TRD/macros/AddTRDefficiency.C++");
      AddTRDefficiency(mgr, "ALL", ci);
      //
      // Resolution
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/TRD/macros/AddTRDresolution.C++");      
      AddTRDresolution(mgr, "ALL", ci);
  }
  //
  // Event Statistics (Jan Fiete)
  //

  if (doEventStat) {
      gROOT->LoadMacro("AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
      AliPhysicsSelection* physSel = physSelTask->GetPhysicsSelection();
      physSel->AddBackgroundIdentification(new AliBackgroundSelection());
  }

  if(doCALO) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/QA/AddTaskCalorimeterQA.C");
      AliAnalysisTaskParticleCorrelation *taskCaloQA = AddTaskCalorimeterQA("ESD", kTRUE, kFALSE);
      taskCaloQA->Dump();
  }

  if(doMUONTrig) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskMTRchamberEfficiency.C");
      AliAnalysisTaskTrigChEff *taskMuonTrig = AddTaskMTRchamberEfficiency();
  }

  if(doMUONEff) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG3/muondep/AddTaskMUONTrackingEfficiency.C");
      AliAnalysisTaskMuonTrackingEff *taskMuonTrackEff = AddTaskMUONTrackingEfficiency();
  }
  
  if (doV0) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskV0QA.C");
      AliAnalysisTaskV0QA *taskv0QA = AddTaskV0QA(kFALSE);
  }


  //
  // Wagons that need MC
  //
  
  
  // Init
  if (!mgr->InitAnalysis()) 
      mgr->PrintStatus();
      mgr->PrintStatus();
  // Run on dataset
  mgr->StartAnalysis("local", chain, 1000);
  timer.Stop();
  timer.Print();
}

