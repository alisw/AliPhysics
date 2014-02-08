// $Id$

AliEmcalParticleMaker* AddTaskEmcalParticleMaker(
  const char *tracksName          = "PicoTracks",
  const char *clustersName        = 0,
  const char *tracksOutName       = "EmcalTracks",
  const char *clustersOutName     = "EmcalClusters",
  const char *taskName            = "AliEmcalParticleMaker"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEmcalParticleMaker", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    ::Error("AddTaskEmcalParticleMaker", "This task requires an input event handler");
    return NULL;
  }

  TString clusName(clustersName);
  if (!clustersName) {
    if (evhand->InheritsFrom("AliESDInputHandler")) {
      ::Info("AddTaskEmcalParticleMaker", "ESD analysis, clustersName = \"CaloClusters\"");
      clusName = "CaloClusters";
    } else {
      ::Info("AddTaskEmcalParticleMaker", "AOD analysis, clustersName = \"caloClusters\"");
      clusName,"caloClusters");
    }
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString tname(taskName);
  tname += "_in_" + tracksName + "_" + clusName + "_out_" + tracksOutName + "_" + clustersOutName;
  AliEmcalParticleMaker *eTask = new AliEmcalParticleMaker(tname);
  eTask->SetTracksName(tracksName);
  eTask->SetClusName(clusName);
  eTask->SetTracksOutName(tracksOutName);
  eTask->SetClusOutName(clustersOutName);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(eTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput(eTask, 0, cinput1);
  
  return eTask;
}
