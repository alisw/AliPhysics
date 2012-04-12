AliEmcalClusTrackMatcherTask* AddTaskEmcalClusTrackMatcher(
						       const char *ntracks        = "Tracks",
						       const char *nclusters      = "CaloClusters",
                                                       const Bool_t DoClusTrack   = 1,
                                                       const Bool_t DoTrackClus   = 0
                                                       )
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalClusTrackMatcher", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalClusTrackMatcher", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliEmcalClusTrackMatcherTask* matcher = new AliEmcalClusTrackMatcherTask();
  matcher->SetTracksName(ntracks);
  matcher->SetClusName(nclusters);
  matcher->SetDoClusTrackMatching(DoClusTrack);
  matcher->SetDoTrackClusMatching(DoTrackClus);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(matcher);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer() ;
  
  mgr->ConnectInput  (matcher, 0,  cinput1 );
  mgr->ConnectOutput (matcher, 0, coutput1 );
  
  return matcher;
  
}
