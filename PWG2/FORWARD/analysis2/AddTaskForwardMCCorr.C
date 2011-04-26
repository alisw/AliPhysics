/**
 * @file   AddTaskForwardMCCorr.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Tue Apr 26 09:56:39 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_forward_scripts_tasks
 */
/** 
 * Add a Forward MC correction generator task to train 
 * 
 * 
 * @return Added task 
 *
 * @ingroup pwg2_forward_mc
 */
AliAnalysisTask*
AddTaskForwardMCCorr()
{
  // --- Creating the manager and handlers ---------------------------
  AliAnalysisManager *mgr  = AliAnalysisManager::GetAnalysisManager();
  if (!mgr->GetMCtruthEventHandler()) { 
    Error("AddTaskCentralMCCorr", 
	  "No MC input handler defined - cannot continue");
    return 0;
  }

  // --- Add our task ------------------------------------------------
  AliForwardMCCorrectionsTask* task = new AliForwardMCCorrectionsTask("fmd");
  mgr->AddTask(task);
  task->GetTrackDensity().SetDebug(false);
  task->GetTrackDensity().SetMaxConsequtiveStrips(3);
  
  // --- create containers for input/output --------------------------
  AliAnalysisDataContainer *sums = 
    mgr->CreateContainer("ForwardSums", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *output = 
    mgr->CreateContainer("ForwardResults", TList::Class(), 
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
