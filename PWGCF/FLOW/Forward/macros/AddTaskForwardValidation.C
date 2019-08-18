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
 //#include "AliForwardTaskValidation.h"

AliAnalysisTaskSE* AddTaskForwardValidation(Bool_t esd, Bool_t mc, Bool_t useEventcuts,TString suffix)
{
  std::cout << "AddTaskForwardValidation" << std::endl;
    AliForwardTaskValidation* task = AliForwardTaskValidation::ConnectTask(suffix, suffix);
    task->fSettings.mc = mc;
    task->fSettings.esd = esd;
    task->fSettings.useEventcuts = useEventcuts;

  return task;
}
/*
 * EOF
 *
 */
