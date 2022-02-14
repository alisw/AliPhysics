/**
 * @file   AddTaskForwardValidation.C
 * @author Freja Thoresen <freja.thoresen@cern.ch>
 *
 * @brief  Add validation forward task to train
 *
 *
 * @ingroup pwglf_forward_scripts_tasks
 */


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
