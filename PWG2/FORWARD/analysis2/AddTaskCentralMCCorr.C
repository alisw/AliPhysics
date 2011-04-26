/**
 * @file   AddTaskCentralMCCorr.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Apr 26 09:55:29 2011
 * 
 * @brief  
 * 
 * @ingroup pwg2_forward_scripts_tasks
 * 
 */

/** 
 * Add a Central MC correction generator task to train 
 * 
 * 
 * @return Added task 
 *
 * @ingroup pwg2_central_mc
 */
AliAnalysisTask*
AddTaskCentralMCCorr()
{
  // --- Creating the manager and handlers ---------------------------
  AliAnalysisManager *mgr  = AliAnalysisManager::GetAnalysisManager();
  if (!mgr->GetMCtruthEventHandler()) { 
    Error("AddTaskCentralMCCorr", 
	  "No MC input handler defined - cannot continue");
    return 0;
  }

  // --- Add our task ------------------------------------------------
  AliCentralMCCorrectionsTask* task2 = new AliCentralMCCorrectionsTask("spd");
  mgr->AddTask(task2);
  task2->SetNPhiBins(40);
  // task2->GetTrackDensity().SetDebug(false);
  
  // --- create containers for input/output --------------------------
  AliAnalysisDataContainer *sums = 
    mgr->CreateContainer("CentralSums", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *output = 
    mgr->CreateContainer("CentralResults", TList::Class(), 
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
