// $Id$

AliEmcalQATask* AddTaskEmcalQA(
  const char *taskname           = "AliEmcalQATask",
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
    ::Error("AddTaskEmcalQA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalQA", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  AliEmcalQATask* qaTask = new AliEmcalQATask(taskname);
  qaTask->SetTracksName(ntracks);
  qaTask->SetClusName(nclusters);
  qaTask->SetJetsName(njets);
  qaTask->SetTrgClusName(ntrgclusters);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(qaTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;

  TString contName(taskname);
  contName += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (qaTask, 0,  cinput1 );
  mgr->ConnectOutput (qaTask, 1, coutput1 );

  return qaTask;
}
