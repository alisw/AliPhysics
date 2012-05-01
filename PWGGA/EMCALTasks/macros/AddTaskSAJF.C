// $Id$

AliAnalysisTaskSAJF* AddTaskSAJF(
  const char *taskname           = "AliAnalysisTaskSAJF",
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  const char *ntrgclusters       = "ClustersL1GAMMAFEE"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskSAJF", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskSAJF", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliAnalysisTaskSAJF* phTask = new AliAnalysisTaskSAJF(taskname);
  phTask->SetTracksName(ntracks);
  phTask->SetClusName(nclusters);
  phTask->SetJetsName(njets);
  phTask->SetTrgClusName(ntrgclusters);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(phTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer()  ;
  TString contname(taskname);
  contname += "_histos";
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(contname.Data(), 
                                                           TList::Class(),AliAnalysisManager::kOutputContainer,
                                                           Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (phTask, 0, cinput);
  mgr->ConnectOutput (phTask, 1, coutput);

  return phTask;
}
