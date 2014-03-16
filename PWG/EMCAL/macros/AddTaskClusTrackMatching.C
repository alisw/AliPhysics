AliAnalysisTaskSE * AddTaskClusTrackMatching(
					     const char*    periodstr          = "LHC11h",
					     const UInt_t   pSel               = AliVEvent::kAny,
					     const char*    inClus             = "EmcCaloClusters",
					     const Double_t trackeff           = 1.0,
					     const Bool_t   doAODTrackProp     = kFALSE,
					     const Bool_t   modifyMatchObjs    = kTRUE,
					     const Bool_t   doHistos           = kFALSE
) {

  // Add task macros for EMCal cluster track matching

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    Error("AddTaskClusTrackMatching","No analysis manager found.");
    return NULL;
  }

  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    Error("AddTaskClusTrackMatching", "This task requires an input event handler");
    return NULL;
  }
  
  TString dType("ESD");
  if (!evhand->InheritsFrom("AliESDInputHandler")) 
    dType = "AOD";

  // Set trackcuts according to period. 
  TString period(periodstr);

  TString inputClus   = TString::Format(inClus);
  TString inputTracks = "AODFilterTracks";

  //----------------------- Filter Tracks -----------------------------------------------------
  const Double_t edist = 440;
  if (dType == "ESD") {
    inputTracks = "ESDFilterTracks";
    TString trackCuts(Form("Hybrid_%s", period.Data()));   
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalEsdTrackFilter.C");
    AliEmcalEsdTrackFilterTask *esdfilter = AddTaskEmcalEsdTrackFilter(inputTracks,trackCuts);
    esdfilter->SetDoPropagation(kTRUE);
    esdfilter->SetDist(edist);
    esdfilter->SelectCollisionCandidates(pSel);
  } else if (dType == "AOD") {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalAodTrackFilter.C");
    AliEmcalAodTrackFilterTask *aodfilter = AddTaskEmcalAodTrackFilter(inputTracks,"tracks",period);
    if (doAODTrackProp) {
      aodfilter->SetDist(edist);
      aodfilter->SetDoPropagation(kTRUE);
    }
    aodfilter->SelectCollisionCandidates(pSel);
  }

  //----------------------- Produce EmcalParticles -----------------------------------------------------
  // Produce objects (AliEmcalParticle) for tracks and clusters 
  // used for cluster-track matching
  TString emctracks = Form("EmcalTracks_%s",inputTracks.Data());
  TString emcclusters = Form("EmcalClusters_%s",inputClus.Data());
  Printf("emctracks: %s  inputTracks: %s",emctracks.Data(),inputTracks.Data());
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalParticleMaker.C");
  AliEmcalParticleMaker *emcalParts = AddTaskEmcalParticleMaker(inputTracks,inputClus,emctracks,emcclusters);
  emcalParts->SelectCollisionCandidates(pSel);

  //----------------------- Cluster-Track matching -----------------------------------------------------
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
  AliEmcalClusTrackMatcherTask *emcalClus =  AddTaskEmcalClusTrackMatcher(emctracks,emcclusters,0.1,modifyMatchObjs,kTRUE);
  emcalClus->SelectCollisionCandidates(pSel);

  Printf("3-- inputTracks: %s emctracks: %s emcclusters: %s",inputTracks.Data(),emctracks.Data(),emcclusters.Data());

  return emcalClus;
}
