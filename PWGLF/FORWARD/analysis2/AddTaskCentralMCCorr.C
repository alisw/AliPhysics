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
 * Add a Central MC correction generator task to train 
 * 
 * 
 * @return Added task 
 *
 * @ingroup pwglf_central_mc
 */
AliAnalysisTask*
AddTaskCentralMCCorr()
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
  AliCentralMCCorrectionsTask* task = new AliCentralMCCorrectionsTask("spd");
  mgr->AddTask(task);
  // This has to match the binning used in the AliAODCentralMult
  // class.  Currently, this is set to 20. 
  task->SetNPhiBins(20);
//  task->SetVertexAxis(40, -20., 20.);
  
  // --- create containers for input/output --------------------------
  AliAnalysisDataContainer *sums = 
    mgr->CreateContainer("CentralCorrSums", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *output = 
    mgr->CreateContainer("CentralCorrResults", TList::Class(), 
			 AliAnalysisManager::kParamContainer, 
			 AliAnalysisManager::GetCommonFileName());

  // --- connect input/output ----------------------------------------
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, sums);
  mgr->ConnectOutput(task, 2, output);

  return task;
}
//
// EOF
// 
