/**
 * @file   AddTaskCopyHeader.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 12:13:43 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_scripts_tasks
 */
/** 
 * Script to add task to copy header from ESD to AOD 
 * 
 * @ingroup pwglf_forward_aod
 */
AliAnalysisTaskSE*
AddTaskCopyHeader(const char* what="")
{
  // --- Load libraries ----------------------------------------------
  gROOT->LoadClass("AliCopyHeaderTask", "libPWGLFforward2");

  // --- Create task -------------------------------------------------
  AliCopyHeaderTask* task = new AliCopyHeaderTask;
  task->SetCopyOptions(what);
  if (!task->Connect()) {
    delete task;
    task = 0;
  }
  return task;
}
//
// EOF
//
