// $Id$

AliEmcalEsdTpcTrackTask* AddTaskEmcalEsdTpcTrack(
  const char *name       = "TpcSpdVertexConstrainedTracks"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalEsdTpcTrack", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalEsdTpcTrack", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/macros/CreateTrackCutsPWGJE.C");

  /* hybrid track cuts*/
  AliESDtrackCuts *cutsp = CreateTrackCutsPWGJE(10001007);       //1000 adds SPD any requirement
  AliESDtrackCuts *hybsp = CreateTrackCutsPWGJE(10041007);       //1004 removes ITSrefit requirement from standard set    
  hybsp->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);

  AliEmcalEsdTpcTrackTask *eTask = new AliEmcalEsdTpcTrackTask();
  eTask->SetTrackCuts(cutsp);
  eTask->SetHybridTrackCuts(hybsp);

  eTask->SetTracksName(name);

  cout << " *** Hybrid track selector task configured *** " << endl;

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(eTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(eTask, 0, cinput1);
  
  return eTask;
}
