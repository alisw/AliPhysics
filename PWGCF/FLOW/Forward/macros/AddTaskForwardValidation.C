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

AliAnalysisTaskSE* AddTaskForwardValidation(TString name,TString suffix)
{
  std::cout << "AddTaskForwardValidation" << std::endl;
    AliForwardTaskValidation* task = AliForwardTaskValidation::ConnectTask(name, suffix);

  return task;
}
/*
 * EOF
 *
 */
