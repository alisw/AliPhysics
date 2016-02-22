AliEmcalClusTrackMatcherTask* AddTaskEmcalClusTrackMatcher(const char *nTracks          = "usedefault",
                                                           const char *nClusters        = "usedefault",
                                                           const Double_t maxDist       = 0.1,
                                                           const Bool_t attachEmcalPart = kFALSE,
                                                           const Bool_t updateClusters  = kTRUE,
                                                           const Bool_t updateTracks    = kTRUE,
                                                           const Bool_t createHisto     = kFALSE)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalClusTrackMatcher", "No analysis manager to connect to.");
    return 0;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalClusTrackMatcher", "This task requires an input event handler");
    return 0;
  }
  
  enum EDataType_t {
    kUnknown,
    kESD,
    kAOD
  };

  EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString trackName(nTracks);
  TString clusName(nClusters);

  if (trackName == "usedefault") {
    if (dataType == kESD) {
      trackName = "Tracks";
    }
    else if (dataType == kAOD) {
      trackName = "tracks";
    }
    else {
      trackName = "";
    }
  }

  if (clusName == "usedefault") {
    if (dataType == kESD) {
      clusName = "CaloClusters";
    }
    else if (dataType == kAOD) {
      clusName = "caloClusters";
    }
    else {
      clusName = "";
    }
  }

  TString name(Form("ClusTrackMatcher_%s_%s", trackName.Data(), clusName.Data()));
  AliEmcalClusTrackMatcherTask* matcher = new AliEmcalClusTrackMatcherTask(name, createHisto);
  AliTrackContainer* trackCont = matcher->AddTrackContainer(trackName);
  if (trackCont) {
    if (trackName == "Tracks" || trackName == "tracks") trackCont->SetFilterHybridTracks(kTRUE);
  }
  matcher->AddClusterContainer(clusName);
  matcher->SetMaxDistance(maxDist);
  matcher->SetUpdateClusters(updateClusters);
  matcher->SetUpdateTracks(updateTracks);

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
