// $Id$

AliAnalysisTaskSE* AddTaskJetPreparation(
  const char*    dataType           = "ESD",
  const char*    periodstr          = "LHC11h",
  const char*    usedTracks         = "PicoTracks",
  const char*    usedMCParticles    = "MCParticles",
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
  const Bool_t   isEmcalTrain       = kFALSE
)
{
  // Add task macros for all jet related helper tasks.

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    Error("AddTaskJetPreparation","No analysis manager found.");
    return 0;
  }
  Bool_t isMC = (mgr->GetMCtruthEventHandler() != NULL);

  // Set trackcuts according to period. Every period used should be definied here
  TString period(periodstr);
  TString clusterColName(usedClusters);
  if ( (period != "LHC11h") && (period!="LHC11a") )
  {
    Error("AddTaskJetPreparation","###################################################");
    Error("AddTaskJetPreparation","Run period in AddTaskJetPreparation.C not recognized! You have to specify it for the used period, if you need jets!");
    Error("AddTaskJetPreparation","###################################################");
    return 0;
  }    

  TString dType(dataType);

  if (dType == "AOD")
    clusterColName = "caloClusters";

  if (makePicoTracks && (dType == "ESD" || dType == "AOD") )
  {
    TString inputTracks = "tracks";

    if (dType == "ESD")
    {
      inputTracks = "HybridTracks";
      TString trackCuts(Form("Hybrid_%s", period.Data()));
      // Hybrid tracks maker for ESD
      gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalEsdTpcTrack.C");
      AliEmcalEsdTpcTrackTask *hybTask = AddTaskEmcalEsdTpcTrack(inputTracks.Data(),trackCuts.Data());
      hybTask->SelectCollisionCandidates(pSel);

      // Track propagator to extend track to the TPC boundaries
      gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalTrackPropagator.C");
      AliEmcalTrackPropagatorTask *propTask = AddTaskEmcalTrackPropagator(inputTracks.Data());
      propTask->SelectCollisionCandidates(pSel);
    }
    // Produce PicoTracks 
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
    AliEmcalPicoTrackMaker *pTrackTask = AddTaskEmcalPicoTrackMaker("PicoTracks", inputTracks.Data(), period.Data());
    pTrackTask->SelectCollisionCandidates(pSel);
  }

  // Produce particles used for hadronic correction
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalParticleMaker.C");
  AliEmcalParticleMaker *emcalParts = AddTaskEmcalParticleMaker(usedTracks,clusterColName.Data(),"EmcalTracks","EmcalClusters");
  emcalParts->SelectCollisionCandidates(pSel);

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
      RequestMemory(hCorr,250*1024);
  }

  // Produce MC particles
  if (isMC && strlen(usedMCParticles)>0)
  {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskMCTrackSelector.C");
    AliEmcalMCTrackSelector *mcPartTask = AddTaskMCTrackSelector(usedMCParticles, kFALSE, kFALSE);
    mcPartTask->SelectCollisionCandidates(pSel);
  }

  // Return one task that represents the jet preparation on LEGO trains
  return emcalParts;
}
