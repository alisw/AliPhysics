// $Id$

AliAnalysisTaskSATR* AddTaskSATR(
  const char* triggerType = "L0",
  AliAnalysisTaskEMCALClusterizeFast *clusterizer = 0, 
  AliAnalysisTaskEMCALClusterizeFast *trgClusterizer = 0)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskSATR", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskSATR", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString caloClusName, trgClusName, taskName;
  if (!strcmp(triggerType, "L0")) {
    caloClusName = "ClustersL0FEE";
    trgClusName = "ClustersL0FOR";
    taskName = "TriggerAnaL0";
  }
  else if (!strcmp(triggerType, "L1GAMMA")) {
    caloClusName = "ClustersL1GAMMAFEE";
    trgClusName = "ClustersL1GAMMAFOR";
    taskName = "TriggerAnaL1GAMMA";
  }
  else if (!strcmp(triggerType, "L1JET")) {
    caloClusName = "ClustersL1JETFEE";
    trgClusName = "ClustersL1JETFOR";
    taskName = "TriggerAnaL1JET";
  }
	
  AliAnalysisTaskSATR* task = new AliAnalysisTaskSATR(taskName.Data());
  task->SetTimeCutOn(kTRUE);
  task->SetL0TimeCut(-20, 20);
  task->SetClusTimeCut(-1, 1);
  task->SetCheckDeadClusters(kFALSE);
  task->SetLoadPed(kFALSE);
  task->SetCaloClustersName(caloClusName.Data());
  task->SetTriggerClustersName(trgClusName.Data());
  task->SetClusterizer(clusterizer);
  task->SetTriggerClusterizer(trgClusterizer);
  task->SetCutL0Amp(-1, 999);
  task->SetCutClusEnergy(-1, 999);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("coutput", 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (task, 0,  cinput1 );
  mgr->ConnectOutput (task, 1, coutput1 );

  return task;
}
