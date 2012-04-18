AliEmcalAodTrackFilterTask* AddTaskEmcalAodTrackFilter(
						       const char *name       = "PicoTracks",
						       const char *inname     = "tracks",
						       const char *runPeriod  = "LHC11h"
                                                       )
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalAodTrackFilter", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalAodTrackFilter", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  // Add aod track filter task.


  AliEmcalAodTrackFilterTask *eTask = new AliEmcalAodTrackFilterTask();
  eTask->SetTracksName(name);
  eTask->SetTracksIn(inname);
  eTask->SetRunPeriod(runPeriod);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(eTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  //AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer() ;
  
  mgr->ConnectInput  (eTask, 0,  cinput1 );
  //mgr->ConnectOutput (eTask, 0, coutput1 );
  
  return eTask;
  
}
