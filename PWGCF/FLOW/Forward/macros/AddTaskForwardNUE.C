/**
 * @file   AddTaskForwardNUA.C
 * @author Freja Thoresen <freja.thoresen@cern.ch>
 *
 * @brief  Add Q-cummulant forward task to train
 *
 *
 */

AliAnalysisTaskSE* AddTaskForwardNUE(Bool_t nua_mode)
{
  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    Fatal("","No analysis manager to connect to.");

  const char* name = Form("ForwardNUE_filled");

  TString resName = "NUE_filled";
  if (nua_mode){
    resName = "NUE_extrapolated";
    name = "ForwardNUE";
  }
  AliForwardNUETask* task = new AliForwardNUETask(name);
  task->nua_mode = nua_mode;
  std::cout << resName << std::endl;

    AliAnalysisDataContainer *coutput_recon =
    mgr->CreateContainer(resName,
     TList::Class(),
     AliAnalysisManager::kOutputContainer,
     mgr->GetCommonFileName());
    task->fSettings.fDataType = task->fSettings.kRECON;
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutput_recon);



  return task;
}
/*
 * EOF
 *
 */
