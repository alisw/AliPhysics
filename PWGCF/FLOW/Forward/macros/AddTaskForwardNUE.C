/**
 * @file   AddTaskForwardNUA.C
 * @author Freja Thoresen <freja.thoresen@cern.ch>
 *
 * @brief  Add Q-cummulant forward task to train
 *
 *
 */

AliAnalysisTaskSE* AddTaskForwardNUE(TString name1)
{
  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    Fatal("","No analysis manager to connect to.");

  const char* name = Form("ForwardNUE_filled");

  TString resName = "NUE_filled";
  AliForwardNUETask* task = new AliForwardNUETask(name);
  std::cout << resName << std::endl;

  AliAnalysisDataContainer* valid = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("event_selection_xchange");
  task->ConnectInput(1,valid);

    AliAnalysisDataContainer *coutput_recon =
    mgr->CreateContainer(resName,
     TList::Class(),
     AliAnalysisManager::kOutputContainer,
     mgr->GetCommonFileName());

    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutput_recon);



  return task;
}
/*
 * EOF
 *
 */
