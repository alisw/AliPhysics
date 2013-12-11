// $Id$

AliAnalysisTaskSE* AddTaskJetPreparation(
  const char*    dataType           = "ESD",
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

  // Set trackcuts according to period. Every period used should be definied here
  TString period(periodstr);
  TString clusterColName(usedClusters);
  TString particleColName(usedMCParticles);
  TString dType(dataType);

  if ((dType == "AOD") && (clusterColName == "CaloClusters"))
    clusterColName = "caloClusters";
  if ((dType == "ESD") && (clusterColName == "caloClusters"))
    clusterColName = "CaloClusters";

  if (makePicoTracks && (dType == "ESD" || dType == "AOD") )
  {
    TString inputTracks = "tracks";

    if (dType == "ESD")
    {
      inputTracks = "HybridTracks";
      TString trackCuts(Form("Hybrid_%s", period.Data()));
      // Hybrid tracks maker for ESD
      gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalEsdTrackFilter.C");
      AliEmcalEsdTrackFilterTask *hybTask = AddTaskEmcalEsdTrackFilter(inputTracks.Data(),trackCuts.Data());
      hybTask->SelectCollisionCandidates(pSel);
      hybTask->SetDoPropagation(kTRUE);
    }
    if(dType == "AOD" && doAODTrackProp) {
      // Track propagator to extend track to the EMCal surface
      gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalTrackPropagatorAOD.C");
      AliEmcalTrackPropagatorTaskAOD *propTask = AddTaskEmcalTrackPropagatorAOD(inputTracks.Data(),440.);
      propTask->SelectCollisionCandidates(pSel);
    }


    // Produce PicoTracks 
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
    AliEmcalPicoTrackMaker *pTrackTask = AddTaskEmcalPicoTrackMaker("PicoTracks", inputTracks.Data(), period.Data());
    pTrackTask->SelectCollisionCandidates(pSel);
    pTrackTask->SetTrackEfficiency(trackeff);
  }

  // Produce particles used for hadronic correction
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalParticleMaker.C");
  AliEmcalParticleMaker *emcalParts = AddTaskEmcalParticleMaker(usedTracks,clusterColName.Data(),"EmcalTracks","EmcalClusters");
  emcalParts->SelectCollisionCandidates(pSel);

  // Trigger maker
  if (makeTrigger) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalTriggerMaker.C");
    AliEmcalTriggerMaker *emcalTriggers = AddTaskEmcalTriggerMaker("EmcalTriggers");
    emcalTriggers->SelectCollisionCandidates(pSel);
  }

  // Relate tracks and clusters
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
  AliEmcalClusTrackMatcherTask *emcalClus =  AddTaskEmcalClusTrackMatcher("EmcalTracks","EmcalClusters",0.1);
  emcalClus->SelectCollisionCandidates(pSel);
  if (isEmcalTrain)
    RequestMemory(emcalClus,100*1024);

  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskHadCorr.C"); 
  AliHadCorrTask *hCorr = AddTaskHadCorr("EmcalTracks","EmcalClusters",outClusName,hadcorr,minPtEt,phiMatch,etaMatch,Eexcl,trackclus,doHistos);
  hCorr->SelectCollisionCandidates(pSel);

  if (isEmcalTrain) {
    if (doHistos)
      RequestMemory(hCorr,500*1024);
  }

  // Produce MC particles
  if(particleColName != "")
  {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskMCTrackSelector.C");
    AliEmcalMCTrackSelector *mcPartTask = AddTaskMCTrackSelector(particleColName.Data(), kFALSE, kFALSE);
    mcPartTask->SelectCollisionCandidates(pSel);
  }

  // Return one task that represents the jet preparation on LEGO trains
  return emcalParts;
}
