AliAnalysisTaskCaloCellsQA* AddTaskCaloCellsPhysQA(Int_t nmods = 10, Int_t det = 0,
                                               TString fname = "CellsQA.root", TString contname = "")
{
  // Task to add EMCAL/PHOS cellsQA/runsQA to your analysis.
  //
  // Usage example for EMCAL:
  //
  //   gSystem->Load("libPWGGAPHOSTasks");
  //   gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/PHOSTasks/CaloCellQA/macros/AddTaskCaloCellsPhysQA.C");
  //   AliAnalysisTaskCaloCellsQA *taskQA = AddTaskCaloCellsPhysQA(10); // 10 supermodules
  //   taskQA->SelectCollisionCandidates(AliVEvent::kMB); // if necessary
  //   // taskQA->SetAvoidPileup(kFALSE); // some customization
  //   // taskQA->GetCaloCellsQA()->ActivateFullAnalysis(); // more histograms, not usually necessary
  //   // Int_t badcells[] = {74,103,917};
  //   // taskQA->SetBadCells(badcells, 3); // reject clusters containing any of these cells
  //
  // Usage example for PHOS:
  //
  //   gSystem->Load("libPWGAPHOSTasks");
  //   gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/PHOSTasks/CaloCellQA/macros/AddTaskCaloCellsPhysQA.C");
  //   AliAnalysisTaskCaloCellsQA *taskQA = AddTaskCaloCellsPhysQA(5, 1);
  //   taskQA->SelectCollisionCandidates(AliVEvent::kMB); // if necessary
  //   taskQA->GetCaloCellsQA()->SetClusterEnergyCuts(0.3,0.1); // increase statistics
  //   // taskQA->SetAvoidPileup(kFALSE); // some customization
  //   // taskQA->GetCaloCellsQA()->ActivateFullAnalysis(); // more histograms, not usually necessary
  //   // Int_t badcells[] = {1234};
  //   // taskQA->SetBadCells(badcells, 1); // reject clusters containing any of these cells
  //
  // nmods -- maximum supermodule number + 1:
  //   use 4 for EMCAL <= 2010;
  //   use 4 for PHOS <=2013 (Run I), use 5 for PHOS >=2015 (PHOS numbers start from 1, not from zero);
  //   use 10 for EMCAL >= 2011;
  // det -- detector, 0/EMCAL, 1/PHOS;
  // fname -- output file name;
  //   if NULL, the output will be written into mgr->GetCommonFileName() + container;
  // contname -- TObjArray container name in the output file;
  //   if not NULL (or fname = NULL), the output will be written into output container and
  //   cannot be later merged for different run numbers;
  //   contname must be unique, if you are going to call AddTaskCaloCellsPhysQA() several times.
  //
  // Note that if fname = NULL and contname = NULL, the output will be written into
  //   file mgr->GetCommonFileName() with container name CellsQAResults.

  // get manager instance

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCaloCellsPhysQA", "No analysis manager to connect to");
    return NULL;
  }

  // check the analysis type using the event handlers connected to the analysis manager
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskCaloCellsPhysQA", "This task requires an input event handler");
    return NULL;
  }

  // Configure analysis
  //===========================================================================

  // detector number
  Int_t det2 = -1;
  if       (det == 0) det2 = AliAnalysisTaskCaloCellsQA::kEMCAL;
  else  if (det == 1) det2 = AliAnalysisTaskCaloCellsQA::kPHOS;
  else
  ::Fatal("AddTaskCaloCellsPhysQA", "Wrong detector provided");

  AliAnalysisTaskCaloCellsQA* task;

  if ((fname.Length() != 0) && (contname.Length() == 0)) task = new AliAnalysisTaskCaloCellsPhysQA("AliAnalysisTaskCaloCellsQA", nmods, det2, fname.Data());
  else                    task = new AliAnalysisTaskCaloCellsPhysQA("AliAnalysisTaskCaloCellsQA", nmods, det2);
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  // container output into particular file
  if ((fname.Length() != 0) && (contname.Length() != 0))
    mgr->ConnectOutput(task, 1, mgr->CreateContainer(contname.Data(),
                             TObjArray::Class(), AliAnalysisManager::kOutputContainer, fname.Data()));

  // container output into common file
  if ((fname.Length() == 0)) {
    if (contname.Length() == 0) contname = "CellsQAResults";
      mgr->ConnectOutput(task, 1, mgr->CreateContainer(contname.Data(),
                               TObjArray::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName()));
  }
  return task;
}