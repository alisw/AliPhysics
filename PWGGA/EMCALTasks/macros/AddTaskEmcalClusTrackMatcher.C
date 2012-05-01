// $Id$

AliEmcalClusTrackMatcherTask* AddTaskEmcalClusTrackMatcher(
  const char *nTracks    = "Tracks",
  const char *nClusters  = "CaloClusters",
  Bool_t doClusTrack     = kTRUE,
  Bool_t doTrackClus     = kFALSE
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
  TString name(Form("ClusTrackMatcher_%s_%s",nTracks,nClusters));
  AliEmcalClusTrackMatcherTask* matcher = new AliEmcalClusTrackMatcherTask(name);//name.Data());
  matcher->SetTracksName(nTracks);
  matcher->SetClusName(nClusters);
  matcher->SetDoClusTrackMatching(doClusTrack);
  matcher->SetDoTrackClusMatching(doTrackClus);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(matcher);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  mgr->ConnectInput  (matcher, 0,  cinput1 );

  return matcher;
}
