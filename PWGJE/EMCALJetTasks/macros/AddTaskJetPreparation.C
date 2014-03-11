// $Id$

AliAnalysisTaskSE* AddTaskJetPreparation(
  const char*    periodstr          = "LHC11h",
  const char*    usedTracks         = "PicoTracks",
  const char*    usedMCParticles    = "MCParticlesSelected",
  const char*    usedClusters       = "CaloClusters",
  const char*    outClusName        = "CaloClustersCorr",
  const Double_t hadcorr            = 2.0,
  const Double_t Eexcl              = 0.00,
  const Double_t phiMatch           = 0.03,
  const Double_t etaMatch           = 0.015,
  const Double_t minPtEt            = 0.15,
  const UInt_t   pSel               = AliVEvent::kAny,
  const Bool_t   trackclus          = kTRUE,
  const Bool_t   doHistos           = kFALSE,
  const Bool_t   makePicoTracks     = kTRUE,
  const Bool_t   makeTrigger        = kTRUE,
  const Bool_t   isEmcalTrain       = kFALSE,
  const Double_t trackeff           = 1.0,
  const Bool_t   doAODTrackProp     = kFALSE
)
{
  // Add task macros for all jet related helper tasks.

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    Error("AddTaskJetPreparation","No analysis manager found.");
    return 0;
  }

  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    Error("AddTaskJetPreparation", "This task requires an input event handler");
    return NULL;
  }

  // Set trackcuts according to period. Every period used should be definied here
  TString period(periodstr);
  TString clusterColName(usedClusters);
  TString particleColName(usedMCParticles);

  TString dType("ESD");
  if (!evhand->InheritsFrom("AliESDInputHandler")) 
    dType = "AOD";
  if ((dType == "AOD") && (clusterColName == "CaloClusters"))
    clusterColName = "caloClusters";
  if ((dType == "ESD") && (clusterColName == "caloClusters"))
    clusterColName = "CaloClusters";

  if (0) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalTrackPropagator.C");
    AliEmcalTrackPropagatorTask *proptask = AddTaskEmcalTrackPropagator();
    proptask->SelectCollisionCandidates(pSel);
  }

  if (makePicoTracks && (dType == "ESD" || dType == "AOD")) {
    // Filter tracks
    TString inputTracks = "AODFilterTracks";
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
    // Produce PicoTracks 
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
    AliEmcalPicoTrackMaker *pTrackTask = AddTaskEmcalPicoTrackMaker("PicoTracks", inputTracks);
    pTrackTask->SetTrackEfficiency(trackeff);
    pTrackTask->SelectCollisionCandidates(pSel);
  }

  // Trigger maker
  if (makeTrigger) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalTriggerMaker.C");
    AliEmcalTriggerMaker *emcalTriggers = AddTaskEmcalTriggerMaker("EmcalTriggers");
    emcalTriggers->SelectCollisionCandidates(pSel);
  }

  TString emctracks(Form("EmcalTracks_%s",usedTracks));
  TString emcclusters(Form("EmcalClusters_%s",clusterColName.Data()));

  // Produce particles used for hadronic correction
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalParticleMaker.C");
  AliEmcalParticleMaker *emcalParts = AddTaskEmcalParticleMaker(usedTracks,clusterColName,emctracks,emcclusters);
  emcalParts->SelectCollisionCandidates(pSel);

  // Relate tracks and clusters
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
  AliEmcalClusTrackMatcherTask *emcalClus =  AddTaskEmcalClusTrackMatcher(emctracks,emcclusters,0.1,kTRUE,doHistos);
  emcalClus->SelectCollisionCandidates(pSel);
  if (isEmcalTrain)
    RequestMemory(emcalClus,100*1024);

  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskHadCorr.C"); 
  AliHadCorrTask *hCorr = AddTaskHadCorr(emctracks,emcclusters,outClusName,hadcorr,
					 minPtEt,phiMatch,etaMatch,Eexcl,trackclus,doHistos);
  hCorr->SelectCollisionCandidates(pSel);
  if (isEmcalTrain) {
    if (doHistos)
      RequestMemory(hCorr,500*1024);
  }

  // Produce MC particles
  if(particleColName != "") {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskMCTrackSelector.C");
    AliEmcalMCTrackSelector *mcPartTask = AddTaskMCTrackSelector(particleColName, kFALSE, kFALSE);
    mcPartTask->SelectCollisionCandidates(pSel);
  }

  // Return one task that represents the jet preparation on LEGO trains
  return emcalParts;
}
