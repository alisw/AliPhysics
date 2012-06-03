// $Id$

AliEmcalParticleMaker* AddTaskEmcalParticleMaker(
  const char *taskName            = "AliEmcalParticleMaker",
  const char *tracksName          = "PicoTracks",
  const char *clustersName        = "CaloClusters",
  const char *tracksOutName       = "EmcalTracks",
  const char *clustersOutName     = "EmcalClusters"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalParticleMaker", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalParticleMaker", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliEmcalParticleMaker *eTask = new AliEmcalParticleMaker(taskName);
  eTask->SetTracksName(tracksName);
  eTask->SetClusName(clustersName);
  eTask->SetTracksOutName(tracksOutName);
  eTask->SetClusOutName(clustersOutName);
  eTask->SetAnaType(AliAnalysisTaskEmcal::kEMCAL);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(eTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput  (eTask, 0,  cinput1 );
  
  return eTask;
}
