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
void
AddTaskCopyHeader()
{
  // --- Load libraries ----------------------------------------------
  gROOT->LoadClass("AliAODForwardMult", "libPWGLFforward2");

  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskCopyHeader", "No analysis manager to connect to.");
    return;
  }   

  // --- Create task -------------------------------------------------
  AliCopyHeaderTask* task = new AliCopyHeaderTask;
  mgr->AddTask(task);
  
  // --- Connect input -----------------------------------------------
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
}
//
// EOF
//
