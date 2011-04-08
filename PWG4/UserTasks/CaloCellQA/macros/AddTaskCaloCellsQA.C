AliAnalysisTaskCaloCellsQA* AddTaskCaloCellsQA(Int_t nmods = 10, Int_t det = 0, char* fname = "CellsQA.root",
                                               Bool_t initGeom = kFALSE, Bool_t kFullAnalysis = kFALSE)
{
  // Task to add EMCAL/PHOS cellsQA/runsQA to your analysis.
  // Do not forget to initialize geometry!
  //
  // Usage example for EMCAL:
  //
  //   gROOT->LoadMacro("$ALICE_ROOT/PWG4/UserTasks/CaloCellQA/macros/AddTaskCaloCellsQA.C");
  //   AliAnalysisTaskCaloCellsQA *taskQA = AddTaskCaloCellsQA(10); // 10 supermodules
  //   taskQA->SelectCollisionCandidates(AliVEvent::kMB); // if necessary
  //   // taskQA->SetAvoidPileup(kFALSE); // some customization
  //
  // Usage example for PHOS:
  //
  //   gROOT->LoadMacro("$ALICE_ROOT/PWG4/UserTasks/CaloCellQA/macros/AddTaskCaloCellsQA.C");
  //   AliAnalysisTaskCaloCellsQA *taskQA = AddTaskCaloCellsQA(4, 1);
  //   taskQA->SelectCollisionCandidates(AliVEvent::kMB); // if necessary
  //   taskQA->GetCaloCellsQA()->SetClusterEnergyCuts(0.3,0.1); // increase statistics
  //
  // fname -- output file name;
  // nmods -- maximum supermodule number + 1:
  //   use 4 for EMCAL <= 2010;
  //   use 4 for PHOS (PHOS numbers start from 1, not from zero);
  //   use 10 for EMCAL >= 2011;
  // det -- detector, 0/EMCAL, 1/PHOS;
  // initGeom -- if true, initialize geometry for you;
  // kFullAnalysis -- if true, initialize the analysis to fill more histograms
  //                  (necessary for the completeness in searching of problematic cells).

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

  AliAnalysisTaskCaloCellsQA* task = new AliAnalysisTaskCaloCellsQA();
  mgr->AddTask(task);

  // initialize geometry
  if (initGeom) {
    if      (det == 0) AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");
    else if (det == 1) AliPHOSGeometry::GetInstance("IHEP");
  }

  // initialize analysis instance
  if (det == 0)// EMCAL
    task->InitCaloCellsQA(fname, nmods, AliAnalysisTaskCaloCellsQA::kEMCAL);
  else if (det == 1)// PHOS
    task->InitCaloCellsQA(fname, nmods, AliAnalysisTaskCaloCellsQA::kPHOS);
  else
    ::Fatal("AddTaskCaloCellsQA", "Wrong detector provided");

  // activate filling of all the histograms
  if (kFullAnalysis)
    task->GetCaloCellsQA()->ActivateFullAnalysis();

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  return task;
}
