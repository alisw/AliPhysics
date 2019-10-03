void runPilot() {
  TStopwatch timer;
  timer.Start();
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/ITS -I$ALICE_PHYSICS -I$ALICE_PHYSICS/TRD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTender");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");
  gSystem->Load("libPWG0selectors");
  gSystem->Load("libPWGPP");
  gSystem->Load("libPWG2");
  gSystem->Load("libPWGmuon");
  gSystem->Load("libPWGmuondep");
  gSystem->Load("libPWG2forward");
  gSystem->Load("libPWG4PartCorrBase");
  gSystem->Load("libPWG4PartCorrDep");
 


  Bool_t doQAsym        = 1;   // (data)
  Bool_t doVZERO        = 1;   // (data)
  Bool_t doVertex       = 1;   // (data)
  Bool_t doSPD          = 1;   // (data/RP)   
  Bool_t doFMD          = 1;   // (data)
  Bool_t doTPC          = 1;   // AliPerformanceTask (data/MCtruth)
  Bool_t doEventStat    = 1;   // (data)
  Bool_t doSDD          = 1;   // (data/RP)
  Bool_t doTRD          = 1;   // TRD 
  Bool_t doITS          = 1;   // ITS
  Bool_t doCALO         = 1;   // Calorimeter
  Bool_t doMUONTrig     = 1;   // MUON trigger
  Bool_t doMUONEff      = 0;   // MUON efficiency  NEEDS geometry
  Bool_t doV0           = 0;   // V0 recosntruction performance NEEDS MCtruth
   
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
  chain->AddFile("/home/morsch/AliRoot/trunk/PWGPP/PilotTrain/data/AliESDs.root");
  

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
      gROOT->LoadMacro("$(ALICE_PHYSICS)/PWGPP/macros/AddTaskVertexESD.C");
      AliAnalysisTaskVertexESD* task3 =  AddTaskVertexESD();
      task3->SelectCollisionCandidates();
  }

  //
  // ITS
  // 
  if (doITS) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/AddTaskPerformanceITS.C");
      AliAnalysisTaskITSTrackingCheck *itsQA = AddTaskPerformanceITS(kFALSE);
  }
  
  //
  // SPD (A. Mastroserio)
  //
  if (doSPD) {
      gROOT->LoadMacro("AddTaskSPDQA.C");
      AliAnalysisTaskSE* task4 = AddTaskSPDQA();
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
      gROOT->LoadMacro("$(ALICE_PHYSICS)/PWGPP/TPC/macros/AddTaskPerformanceTPCQA.C");
      AliPerformanceTask *tpcQA = AddTaskPerformanceTPCQA(kFALSE, kTRUE);
  }
  
  //
  // TRD (Alex Bercuci, M. Fasel) 
  //
  if(doTRD) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/AddTrainPerformanceTRD.C");
    AddTrainPerformanceTRD("ALL");
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
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG4/macros/QA/AddTaskCalorimeterQA.C");
      AliAnalysisTaskParticleCorrelation *taskCaloQA = AddTaskCalorimeterQA("ESD", kTRUE, kFALSE);
      taskCaloQA->SetDebugLevel(0);
  }

  if(doMUONTrig) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/AddTaskMTRchamberEfficiency.C");
      AliAnalysisTaskTrigChEff *taskMuonTrig = AddTaskMTRchamberEfficiency();
  }

  if(doMUONEff) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG3/muondep/AddTaskMUONTrackingEfficiency.C");
      AliAnalysisTaskMuonTrackingEff *taskMuonTrackEff = AddTaskMUONTrackingEfficiency();
  }
  
  if (doV0) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/AddTaskV0QA.C");
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
  mgr->StartAnalysis("local", chain);
  timer.Stop();
  timer.Print();
}

