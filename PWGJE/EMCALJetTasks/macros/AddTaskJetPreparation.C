// $Id$

AliAnalysisTaskSE* AddTaskJetPreparation(
  const char*    periodstr          = "LHC11h",
  const char*    pTracksName        = "PicoTracks",
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
  const Bool_t   doAODTrackProp     = kFALSE,
  const Bool_t   modifyMatchObjs    = kTRUE
)
{
  // Add task macros for all jet related helper tasks.

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    Error("AddTaskJetPreparation","No analysis manager found.");
    return NULL;
  }

  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    Error("AddTaskJetPreparation", "This task requires an input event handler");
    return NULL;
  }

  // Set trackcuts according to period. Every period used should be defined here
  TString period(periodstr);
  TString clusterColName(usedClusters);
  TString particleColName(usedMCParticles);
  TString picoTracksName(pTracksName);

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

  //----------------------- Trigger Maker -----------------------------------------------------
  if (makeTrigger) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalTriggerMaker.C");
    AliEmcalTriggerMaker *emcalTriggers = AddTaskEmcalTriggerMaker("EmcalTriggers");
    emcalTriggers->SelectCollisionCandidates(pSel);
  }

  //----------------------- Track Matching tasks -----------------------------------------------------
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskMatchingChain.C");
  AliEmcalClusTrackMatcherTask *emcalClus =  AddTaskMatchingChain(periodstr,pSel,
								  clusterColName,
								  trackeff,doAODTrackProp,
								  0.1,modifyMatchObjs,doHistos);
  
  //hard coded names of AliEmcalParticle strings to coincide with AddTaskClusTrackMatching
  TString inputTracks = "AODFilterTracks";
  if (dType == "ESD") inputTracks = "ESDFilterTracks";
  TString emctracks = Form("EmcalTracks_%s",inputTracks.Data());
  TString emcclusters = Form("EmcalClusters_%s",clusterColName.Data());
  Printf("1-- inputTracks: %s, emcclusters: %s, emctracks: %s",inputTracks.Data(),emcclusters.Data(),emctracks.Data());
  if(makePicoTracks) {
    //----------------------- Produce PicoTracks -----------------------------------------------------
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
    AliEmcalPicoTrackMaker *pTrackTask = AddTaskEmcalPicoTrackMaker(picoTracksName, inputTracks);
    //    pTrackTask->SetTrackEfficiency(trackeff); //now done in Esd/AodFilter
    pTrackTask->SelectCollisionCandidates(pSel);
  }

  //----------------------- Hadronic Correction -----------------------------------------------------
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
  return hCorr;
}
