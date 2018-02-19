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
  Int_t mode = kRECON;

  std::cout << "AddTaskForwardNUA" << std::endl;

  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    Fatal("","No analysis manager to connect to.");

  const char* name = Form("ForwardNUA");
  AliForwardNUATask* task = new AliForwardNUATask(name);

  TString resName = "awesome";


  
  else {
    task->fSettings.fNRefEtaBins = 1; // eller skal det være et andet antal?
  }

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
