AliAnalysisTaskCaloCellsQA* AddTaskCaloCellsQA(Int_t nmods = 10, Int_t det = 0,
                                               char* fname = "CellsQA.root", char* contname = NULL)
{
  // Task to add EMCAL/PHOS cellsQA/runsQA to your analysis.
  //
  // Usage example for EMCAL:
  //
  //   gSystem->Load("libPWG4UserTasks.so");
  //   gROOT->LoadMacro("$ALICE_ROOT/PWG4/UserTasks/CaloCellQA/macros/AddTaskCaloCellsQA.C");
  //   AliAnalysisTaskCaloCellsQA *taskQA = AddTaskCaloCellsQA(10); // 10 supermodules
  //   taskQA->SelectCollisionCandidates(AliVEvent::kMB); // if necessary
  //   // taskQA->SetAvoidPileup(kFALSE); // some customization
  //   // taskQA->GetCaloCellsQA()->ActivateFullAnalysis(); // more histograms, not usually necessary
  //   // Int_t badcells[] = {74,103,917};
  //   // taskQA->SetBadCells(badcells, 3); // reject clusters containing any of these cells
  //
  // Usage example for PHOS:
  //
  //   gSystem->Load("libPWG4UserTasks.so");
  //   gROOT->LoadMacro("$ALICE_ROOT/PWG4/UserTasks/CaloCellQA/macros/AddTaskCaloCellsQA.C");
  //   AliAnalysisTaskCaloCellsQA *taskQA = AddTaskCaloCellsQA(4, 1);
  //   taskQA->SelectCollisionCandidates(AliVEvent::kMB); // if necessary
  //   taskQA->GetCaloCellsQA()->SetClusterEnergyCuts(0.3,0.1); // increase statistics
  //   // taskQA->SetAvoidPileup(kFALSE); // some customization
  //   // taskQA->GetCaloCellsQA()->ActivateFullAnalysis(); // more histograms, not usually necessary
  //   // Int_t badcells[] = {1234};
  //   // taskQA->SetBadCells(badcells, 1); // reject clusters containing any of these cells
  //
  // nmods -- maximum supermodule number + 1:
  //   use 4 for EMCAL <= 2010;
  //   use 4 for PHOS (PHOS numbers start from 1, not from zero);
  //   use 10 for EMCAL >= 2011;
  // det -- detector, 0/EMCAL, 1/PHOS;
  // fname -- output file name;
  //   if NULL, the output will be written into mgr->GetCommonFileName() + container;
  // contname -- TObjArray container name in the output file;
  //   if not NULL (or fname = NULL), the output will be written into output container and
  //   cannot be later merged for different run numbers;
  //   contname must be unique, if you are going to call AddTaskCaloCellsQA() several times.
  //
  // Note that if fname = NULL and contname = NULL, the output will be written into
  //   file mgr->GetCommonFileName() with container name CellsQAResults.

  // get manager instance
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCaloCellsQA", "No analysis manager to connect to");
    return NULL;
  }

  // check the analysis type using the event handlers connected to the analysis manager
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskCaloCellsQA", "This task requires an input event handler");
    return NULL;
  }

  // Configure analysis
  //===========================================================================

  // detector number
  Int_t det2 = -1;
  if       (det == 0) det2 = AliAnalysisTaskCaloCellsQA::kEMCAL;
  else  if (det == 1) det2 = AliAnalysisTaskCaloCellsQA::kPHOS;
  else
  ::Fatal("AddTaskCaloCellsQA", "Wrong detector provided");

  AliAnalysisTaskCaloCellsQA* task;

  if (fname && !contname) task = new AliAnalysisTaskCaloCellsQA("AliAnalysisTaskCaloCellsQA", nmods, det2, fname);
  else                    task = new AliAnalysisTaskCaloCellsQA("AliAnalysisTaskCaloCellsQA", nmods, det2);
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  // container output into particular file
  if (fname && contname)
    mgr->ConnectOutput(task, 1, mgr->CreateContainer(contname,
                       TObjArray::Class(), AliAnalysisManager::kOutputContainer, fname));

  // container output into common file
  if (!fname) {
    if (!contname) contname = "CellsQAResults";
    mgr->ConnectOutput(task, 1, mgr->CreateContainer(contname,
                       TObjArray::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName()));
  }

  return task;
}
