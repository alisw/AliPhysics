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

AliAnalysisTaskSE* AddTaskForwardValidation("name")
{
  std::cout << "AddTaskForwardValidation" << std::endl;

  return AliForwardTaskValidation::ConnectTask("ftvalid", true);
  //std::cout << validation_task << std::endl;
  //return validation_task;
}
/*
 * EOF
 *
 */
