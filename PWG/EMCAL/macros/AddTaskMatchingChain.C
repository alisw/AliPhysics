AliAnalysisTaskSE * AddTaskMatchingChain(
					 const char*    periodstr          = "LHC11h",
					 const UInt_t   pSel               = AliVEvent::kAny,
					 const char*    inClus             = "EmcCaloClusters",
					 const Double_t trackeff           = 1.0,
					 const Bool_t   doAODTrackProp     = kTRUE,
					 const Double_t maxMatchR          = 0.1,
					 const Bool_t   modifyMatchObjs    = kTRUE,
					 const Bool_t   doHistos           = kFALSE,
                                         const Int_t    nCentBins          = 4
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
#ifdef __CLING__
    std::stringstream esdfilteradd;
    esdfilteradd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWG/EMCAL/macros/AddTaskEmcalEsdTrackFilter.C("
        << "\"" << inputTracks << "\", \"" << trackCuts << "\")";
    std::string esdfilteraddstring = esdfilteradd.str();
    AliEmcalEsdTrackFilterTask *esdfilter = (AliEmcalEsdTrackFilterTask *)gROOT->ProcessLine(esdfilteraddstring.c_str());
#else
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalEsdTrackFilter.C");
    AliEmcalEsdTrackFilterTask *esdfilter = AddTaskEmcalEsdTrackFilter(inputTracks,trackCuts);
#endif
    esdfilter->SetDoPropagation(kTRUE);
    esdfilter->SetDist(edist);
    esdfilter->SelectCollisionCandidates(pSel);
    esdfilter->SetTrackEfficiency(trackeff);
  } else if (dType == "AOD") {
#ifdef __CLING__
    std::stringstream aodfilteradd;
    aodfilteradd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWG/EMCAL/macros/AddTaskEmcalAodTrackFilter.C("
        << "\"" << inputTracks << "\", \"tracks\", " << period << "\")";
    std::string aodfilteraddstring = aodfilteradd.str();
    AliEmcalAodTrackFilterTask *aodfilter = (AliEmcalAodTrackFilterTask *)gROOT->ProcessLine(aodfilteraddstring.c_str());
#else
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalAodTrackFilter.C");
    AliEmcalAodTrackFilterTask *aodfilter = AddTaskEmcalAodTrackFilter(inputTracks,"tracks",period);
#endif
    if (doAODTrackProp) {
      aodfilter->SetDist(edist);
      aodfilter->SetAttemptPropMatch(kTRUE);
    }
    aodfilter->SelectCollisionCandidates(pSel);
    aodfilter->SetTrackEfficiency(trackeff);
  }

  //----------------------- Produce EmcalParticles -----------------------------------------------------
  // Produce objects (AliEmcalParticle) for tracks and clusters 
  // used for cluster-track matching
  TString emctracks = Form("EmcalTracks_%s",inputTracks.Data());
  TString emcclusters = Form("EmcalClusters_%s",inputClus.Data());
  Printf("emctracks: %s  inputTracks: %s",emctracks.Data(),inputTracks.Data());
#ifdef __CLING__
  std::stringstream particlemakeradd;
  particlemakeradd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWG/EMCAL/macros/AddTaskEmcalParticleMaker.C("
      << "\"" << inputTracks << "\", \"" << inputClus << "\", \"" << emctracks << "\", \"" << emcclusters << "\")";
  std::string particlemakeraddstring = particlemakeradd.str();
  AliEmcalParticleMaker *emcalParts = (AliEmcalParticleMaker *)gROOT->ProcessLine(particlemakeraddstring.c_str());
#else
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalParticleMaker.C");
  AliEmcalParticleMaker *emcalParts = AddTaskEmcalParticleMaker(inputTracks,inputClus,emctracks,emcclusters);
#endif
  emcalParts->SelectCollisionCandidates(pSel);
  emcalParts->SetNCentBins(nCentBins);
  //----------------------- Cluster-Track matching -----------------------------------------------------
#ifdef __CLING__
  std::stringstream clustertrackmatcheradd;
  clustertrackmatcheradd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C("
      << "\"" << emctracks << "\", \"" << emcclusters << "\", " << maxMatchR << ", " << (modifyMatchObjs ? "kTRUE" : "kFALSE")
      << ", " << (doHistos ? "kTRUE" : "kFALSE") << ")";
  std::string clustertrackmatcheraddstring = clustertrackmatcheradd.str();
  AliEmcalClusTrackMatcherTask *emcalClus = (AliEmcalClusTrackMatcherTask *)gROOT->ProcessLine(clustertrackmatcheraddstring.c_str());
#else
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
  AliEmcalClusTrackMatcherTask *emcalClus =  AddTaskEmcalClusTrackMatcher(emctracks,emcclusters,maxMatchR,modifyMatchObjs,doHistos);
#endif
  emcalClus->SelectCollisionCandidates(pSel);
  emcalClus->SetNCentBins(nCentBins);

  Printf("3-- inputTracks: %s emctracks: %s emcclusters: %s",inputTracks.Data(),emctracks.Data(),emcclusters.Data());

  return emcalClus;
}
