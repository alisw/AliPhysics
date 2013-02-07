// $Id: AddTaskSAQA.C 60163 2013-01-03 09:37:04Z loizides $

AliAnalysisTaskCLQA* AddTaskCLQA(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  const char *taskname           = "AliAnalysisTaskSAQA"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskSAQA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskSAQA", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  
  TString name(taskname);
  if (strcmp(ntracks,"")) {
    name += "_";
    name += ntracks;
  }
  if (strcmp(nclusters,"")) {
    name += "_";
    name += nclusters;
  }
  if (strcmp(njets,"")) {
    name += "_";
    name += njets;
  }

  AliAnalysisTaskCLQA* qaTask = new AliAnalysisTaskCLQA(name);
  qaTask->SetTracksName(ntracks);
  qaTask->SetClusName(nclusters);
  qaTask->SetJetsName(njets);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(qaTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer()  ;

  TString contName(name);
  contName += "_histos";
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(contName.Data(), 
							   TList::Class(),
							   AliAnalysisManager::kOutputContainer,
							   AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (qaTask, 0,  cinput );
  mgr->ConnectOutput (qaTask, 1, coutput );

  return qaTask;
}
