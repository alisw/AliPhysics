// $Id$

//#define EMCALTRAIN

void AddTaskJetPreparation(
  const char*    dataType           = "ESD",
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
  const Bool_t   doHistos           = kFALSE
)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    Error("AddTaskJetPreparation","No analysis manager found.");
    return 0;
  }
  Bool_t isMC = (mgr->GetMCtruthEventHandler() != NULL);

  // Set trackcuts according to period. Every period used should be definied here
  TString period("");
  TString clusterColName(usedClusters);
  TString eper(gSystem->Getenv("ETRAIN_PERIOD"));
  if (eper.BeginsWith("lhc12g") || eper.BeginsWith("lhc11h") || eper.BeginsWith("lhc12a15a") || eper.BeginsWith("LHC12g") || eper.BeginsWith("LHC12a15a") || eper.BeginsWith("lhc13b"))
    period = "LHC11h";
  else if (eper.BeginsWith("lhc11a"))
    period = "LHC11a";
  else
  {
    Error("AddTaskJetPreparation","###################################################");
    Error("AddTaskJetPreparation","Run period in AddTaskJetPreparation.C not recognized! You have to specify it for the used period, if you need jets!");
    return 0;
  }    
 
  Bool_t makePicoTracks = kTRUE;
  if ((eper == "lhc10hs") || (eper == "lhc11hs")) {
    makePicoTracks = kFALSE;
  }

  if((strcmp(dataType,"AOD") == 0) && (clusterColName == "CaloClusters"))
    clusterColName = "caloClusters";

  if( makePicoTracks && ((strcmp(dataType,"ESD") == 0) || (strcmp(dataType,"AOD") == 0)) )
  {
    TString inputTracks = "tracks";

    if(strcmp(dataType,"ESD") == 0)
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
    // PicoTracks maker to produce pico tracks
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
    AliEmcalPicoTrackMaker *pTrackTask = AddTaskEmcalPicoTrackMaker("PicoTracks", inputTracks.Data(), period.Data());
    pTrackTask->SelectCollisionCandidates(pSel);
  }

  // Make particles used for hadronic correction
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalParticleMaker.C");
  AliEmcalParticleMaker *emcalParts = AddTaskEmcalParticleMaker(usedTracks,clusterColName.Data(),"EmcalTracks","EmcalClusters");
  emcalParts->SelectCollisionCandidates(pSel);

  // Relate tracks and clusters
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
  AliEmcalClusTrackMatcherTask *emcalClus =  AddTaskEmcalClusTrackMatcher("EmcalTracks","EmcalClusters",0.1);
  emcalClus->SelectCollisionCandidates(pSel);
  #ifdef EMCALTRAIN
    RequestMemory(emcalClus,100*1024);
  #endif
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskHadCorr.C"); 
  AliHadCorrTask *hCorr = AddTaskHadCorr("EmcalTracks","EmcalClusters",outClusName,hadcorr,minPtEt,phiMatch,etaMatch,Eexcl,doHistos);
  hCorr->SelectCollisionCandidates(pSel);
  #ifdef EMCALTRAIN
    if (doHistos)
      RequestMemory(hCorr,250*1024);
  #endif

  if (isMC && strlen(usedMCParticles)>0)
  {
    // Make MC particles
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskMCTrackSelector.C");
    AliEmcalMCTrackSelector *mcPartTask = AddTaskMCTrackSelector(usedMCParticles, kFALSE, kFALSE);
    mcPartTask->SelectCollisionCandidates(AliVEvent::kAny);
  }
}
