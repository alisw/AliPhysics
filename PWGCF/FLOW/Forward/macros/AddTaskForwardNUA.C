/**
 * @file   FTAddMyTask.C
 * @author Freja Thoresen <freja.thoresen@cern.ch>
 * 
 * @brief  Add Q-cummulant forward task to train 
 * 
 * 
 * @ingroup pwglf_forward_scripts_tasks
 */
/** 
 * @defgroup pwglf_forward_flow Flow 
 *
 * Code to deal with flow 
 *
 * @ingroup pwglf_forward_topical
 */
/** 
 * Add Flow task to train 
 * 
 * @ingroup pwglf_forward_flow
 */
AliAnalysisTaskSE* AddTaskForwardNUA()
{
  Bool_t etagap = true;
  Int_t mode = kRECON;

  std::cout << "AddTaskForwardNUA" << std::endl;

  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    Fatal("","No analysis manager to connect to.");

  const char* name = Form("ForwardNUA");
  AliForwardNUATask* task = new AliForwardNUATask(name);

  TString resName = "awesome";

  if (etagap){
    // if etagap otherwise comment out, and it will be standard
    task->fSettings.fFlowFlags = task->fSettings.kEtaGap;
    task->fSettings.fNRefEtaBins = 1;
  }
  else {
    task->fSettings.fNRefEtaBins = 1; // eller skal det være et andet antal?
  }
  task->fSettings.fUseFMD = true;
  // task->fSettings.fUseV0 = true;
  // V0 has only 8 segments in phi
  if (task->fSettings.fUseFMD) {
    task->fSettings.fNPhiBins = 20;
  } else if (task->fSettings.fUseV0) {
    task->fSettings.fNPhiBins = 8;
  }
  
  task->fSettings.fUseSPDtracklets = true;

  task->fSettings.fZVtxAcceptanceLowEdge = -10;
  task->fSettings.fZVtxAcceptanceUpEdge = 10;
  task->fSettings.fNZvtxBins = 20;

  // Remember to disable multselection framework if necessary!
  task->fSettings.fMultEstimator = "V0M";// RefMult08; // "V0M" // "SPDTracklets";
  // task->fSettings.fMultEstimator = task->fSettings.fMultEstimatorValidTracks;


  if (mode == kRECON) {
    AliAnalysisDataContainer *coutput_recon =
    mgr->CreateContainer(resName,
     TList::Class(),
     AliAnalysisManager::kOutputContainer,
     mgr->GetCommonFileName());
    task->fSettings.fDataType = task->fSettings.kRECON;
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutput_recon);
  }
  else if (mode == kTRUTH) {
    AliAnalysisDataContainer *coutput_truth =
    mgr->CreateContainer(resName,
     TList::Class(),
     AliAnalysisManager::kOutputContainer,
     mgr->GetCommonFileName());

    task->fSettings.fDataType = task->fSettings.kMCTRUTH;
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutput_truth);
  }
  else {
    ::Error("AddTaskForwardNUA", "Invalid mode specified");
  }


  return task;
}
/*
 * EOF
 *
 */
