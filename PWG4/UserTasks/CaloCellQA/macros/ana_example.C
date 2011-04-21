void ana_example()
{
  // Example macro to run QA

  // load relevant library
  gSystem->Load("libPWG4UserTasks.so");

  // change next line to a working code
  TChain* chain = NULL; //CreateChain("wn.xml");
  if (!chain) {
    fprintf(stderr, "FATAL: chain is NULL\n");
    abort();
  }

  AliAnalysisManager *mgr = new AliAnalysisManager("Manager");

  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  // event selection task
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();

  // EMCAL
  Int_t nbadEMCAL = 44;
  Int_t badcellsEMCAL[44] = {
    74,152,495,871,917,1059,1263,1275,1276,1288,1376,1384,1519,1712,1967,2026,2112,
    2114,2115,2116,2117,2120,2123,2298,2671,2768,2769,2770,2771,2773,2774,2776,
    2777,2778,2779,2780,2783,3544,3567,
    103,1382,1961,2047,2540
  };

  // with PS + no pileup
  AliAnalysisTaskCaloCellsQA *task1 = AddTaskCaloCellsQA(4, 0, "CellsQAEMCAL.root");
  task1->SelectCollisionCandidates(AliVEvent::kMB);
  task1->SetBadCells(badcellsEMCAL, nbadEMCAL);

  // with PS + with pileup
  AliAnalysisTaskCaloCellsQA *task2 = AddTaskCaloCellsQA(4, 0, "CellsQAEMCAL2.root");
  task2->SelectCollisionCandidates(AliVEvent::kMB);
  task2->SetBadCells(badcellsEMCAL, nbadEMCAL);
  task2->SetAvoidPileup(kFALSE);

  // no PS + with pileup
  AliAnalysisTaskCaloCellsQA *task3 = AddTaskCaloCellsQA(4, 0, "CellsQAEMCAL3.root");
  task3->SetBadCells(badcellsEMCAL, nbadEMCAL);
  task3->SetAvoidPileup(kFALSE);

  // PHOS
  Int_t nbadPHOS = 1;
  Int_t badcellsPHOS[1] = {9938};

  // with PS + no pileup
  AliAnalysisTaskCaloCellsQA *task4 = AddTaskCaloCellsQA(4, 1, "CellsQAPHOS.root");
  task4->SelectCollisionCandidates(AliVEvent::kMB);
  task4->GetCaloCellsQA()->SetClusterEnergyCuts(0.3,0.1); // increase statistics
  task4->SetBadCells(badcellsPHOS, nbadPHOS);

  // with PS + with pileup
  AliAnalysisTaskCaloCellsQA *task5 = AddTaskCaloCellsQA(4, 1, "CellsQAPHOS2.root");
  task5->SelectCollisionCandidates(AliVEvent::kMB);
  task5->GetCaloCellsQA()->SetClusterEnergyCuts(0.3,0.1);
  task5->SetBadCells(badcellsPHOS, nbadPHOS);
  task5->SetAvoidPileup(kFALSE);

  // no PS + with pileup
  AliAnalysisTaskCaloCellsQA *task6 = AddTaskCaloCellsQA(4, 1, "CellsQAPHOS3.root");
  task6->GetCaloCellsQA()->SetClusterEnergyCuts(0.3,0.1);
  task6->SetBadCells(badcellsPHOS, nbadPHOS);
  task6->SetAvoidPileup(kFALSE);


  if (!mgr->InitAnalysis()) abort();
  mgr->StartAnalysis("local", chain);
}
