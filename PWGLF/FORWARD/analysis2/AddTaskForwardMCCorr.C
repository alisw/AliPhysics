/**
 * @file   AddTaskForwardMCCorr.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Apr 26 09:56:39 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_scripts_tasks
 */
/** 
 * Add a Forward MC correction generator task to train 
 * 
 * 
 * @return Added task 
 *
 * @ingroup pwglf_forward_mc
 */
AliAnalysisTask*
AddTaskForwardMCCorr()
{
  // --- Load libraries ----------------------------------------------
  gROOT->LoadClass("AliAODForwardMult", "libPWGLFforward2");

  // --- Creating the manager and handlers ---------------------------
  AliAnalysisManager *mgr  = AliAnalysisManager::GetAnalysisManager();
  if (!mgr->GetMCtruthEventHandler()) { 
    Error("AddTaskCentralMCCorr", 
	  "No MC input handler defined - cannot continue");
    return 0;
  }

  // --- Add our task ------------------------------------------------
  AliForwardMCCorrectionsTask* task = 
    new AliForwardMCCorrectionsTask("ForwardCorr");
  task->GetTrackDensity().SetDebug(false);
  AliFMDMCTrackDensity& dn = 
    static_cast<AliFMDMCTrackDensity&>(task->GetTrackDensity());
  dn.SetMaxConsequtiveStrips(3);
  //  task->SetVertexAxis(40, -20., 20.);
  
  // --- connect input/output ----------------------------------------
  task->Connect(0, 0);

  return task;
}
//
// EOF
// 
