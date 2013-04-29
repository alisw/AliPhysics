// $Id$

AliEmcalClusTrackMatcherTask* AddTaskEmcalClusTrackMatcher(
  const char *nTracks    = "EmcalTracks",
  const char *nClusters  = "EmcalClusters",
  Double_t maxDist       = 0.1,
  Bool_t   createHisto   = kFALSE
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
  AliEmcalClusTrackMatcherTask* matcher = new AliEmcalClusTrackMatcherTask(name, createHisto);
  matcher->SetTracksName(nTracks);
  matcher->SetClusName(nClusters);
  matcher->SetMaxDistance(maxDist);
  matcher->SetAnaType(AliAnalysisTaskEmcal::kEMCAL);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(matcher);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  mgr->ConnectInput  (matcher, 0,  cinput1 );

  if (createHisto) {
    TString contname(name);
    contname += "_histos";
    AliAnalysisDataContainer *coutput = mgr->CreateContainer(contname,
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     AliAnalysisManager::GetCommonFileName());
    mgr->ConnectOutput(matcher,1,coutput);
  }

  return matcher;
}
