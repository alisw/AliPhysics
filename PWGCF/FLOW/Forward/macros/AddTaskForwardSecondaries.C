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
 #include "AliForwardSettings.h"
 #include "AliAnalysisDataContainer.h"
 #include "AliAnalysisDataSlot.h"

AliAnalysisTaskSE* AddTaskForwardSecondaries()
{
  std::cout << "AddTaskForwardSecondaries" << std::endl;


  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    Fatal("","No analysis manager to connect to.");

  const char* name = Form("ForwardFlowQC");
  AliForwardSecondariesTask* task = new AliForwardSecondariesTask(name);
  TString resName = "Secondaries";


  task->fSettings.fileName = resName;
  mgr->AddTask(task);

  AliAnalysisDataContainer *coutput_recon =
  mgr->CreateContainer(resName,
   AliForwardFlowResultStorage::Class(), //TList::Class(),
   AliAnalysisManager::kOutputContainer,
   mgr->GetCommonFileName());
  //task->fSettings.fDataType = task->fSettings.kRECON;
  mgr->ConnectOutput(task, 1, coutput_recon);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  AliAnalysisDataContainer* valid = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("event_selection_xchange");
  task->ConnectInput(1,valid);


  return task;
}
/*
 * EOF
 *
 */
