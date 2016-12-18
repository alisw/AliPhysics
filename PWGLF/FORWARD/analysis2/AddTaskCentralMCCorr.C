/**
 * @file   AddTaskCentralMCCorr.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Apr 26 09:55:29 2011
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_scripts_tasks
 * 
 */

/** 
 * Add a Central MC correction generator task to train.  
 * This task generates corrections to be stored in OADB file 
 * @c spd_corrections.root 
 * 
 * @return Added task 
 *
 * @ingroup pwglf_central_mc
 */
AliAnalysisTask*
AddTaskCentralMCCorr(Bool_t satellite=false, Bool_t eff=false)
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
  AliCentralMCCorrectionsTask* task = 
    new AliCentralMCCorrectionsTask("CentralCorr");
  // This has to match the binning used in the AliAODCentralMult
  // class.  Currently, this is set to 20. 
  task->SetSatellite(satellite);
  task->SetNPhiBins(20);
  task->SetEffectiveCorrection(eff);
 //  task->SetVertexAxis(40, -20., 20.);
  
  // --- create containers for input/output --------------------------
  task->Connect(0,0);

  return task;
}
//
// EOF
// 
