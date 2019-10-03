// $Id$

AliAnalysisTaskSE *AddTaskJetPreparationHF(
                                           const char*    periodstr          = "LHC11h",
                                           const char*    pTracksName        = "PicoTracks",
                                           Bool_t         isMC               = kFALSE,
                                           const char*    usedMCParticles    = "MCParticlesSelected",
                                           const char*    usedClusters       = "CaloClusters",
                                           const char*    outClusName        = "CaloClustersCorr",
                                           Double_t hadcorr                  = 2.0,
                                           Double_t Eexcl                    = 0.00,
                                           Double_t phiMatch                 = 0.03,
                                           Double_t etaMatch                 = 0.015,
                                           Double_t minPtEt                  = 0.15,
                                           UInt_t   pSel                     = AliVEvent::kAny,
                                           Bool_t   trackclus                = kTRUE,
                                           Bool_t   doHistos                 = kFALSE,
                                           Bool_t   makePicoTracks           = kTRUE,
                                           Bool_t   makeTrigger              = kTRUE,
                                           Bool_t   isEmcalTrain             = kFALSE,
                                           Double_t trackeff                 = 1.0,
                                           Bool_t   doAODTrackProp           = kTRUE,
                                           Bool_t   modifyMatchObjs          = kTRUE,
                                           Bool_t   useOldBitConfig          = kFALSE,
                                           Bool_t   doTriggerQA              = kFALSE,
                                           Int_t    nCentBins                = 4
                                           )
{
  // Add task macros for all jet related helper tasks.

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskJetPreparation", "No analysis manager found.");
    return NULL;
  }

  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    Error("AddTaskJetPreparation", "This task requires an input event handler");
    return NULL;
  }

  Bool_t isESD = kTRUE;
  if (!evhand->InheritsFrom("AliESDInputHandler"))
    isESD = kFALSE; // is AOD    dType = "AOD";

//  // Set trackcuts according to period. Every period used should be defined here
//  TString period(periodstr);
//  TString clusterColName(usedClusters);
//  TString particleColName(usedMCParticles);
//  TString picoTracksName(pTracksName);

  //TString dType("ESD");
  //  if ((dType == "AOD") && (clusterColName == "CaloClusters"))
//    clusterColName = "caloClusters";
//  if ((dType == "ESD") && (clusterColName == "caloClusters"))
//    clusterColName = "CaloClusters";

  //----------------------- Trigger Maker -----------------------------------------------------
  if (makeTrigger) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalTriggerMaker.C");
    /*
     * Parameters (with default values):
     *   triggersOutName      (const char *)       = "EmcalTriggers",
     *   triggerSetupOutName  (const char *)       = "EmcalTriggerSetup",
     *   cellsName            (const char *)       = 0,
     *   triggersName         (const char *)       = 0,
     *   taskName             (const char *)       = "AliEmcalTriggerMaker",
     *   jetLowA              (int)                = 0,
     *   jetLowB              (int)                = 0,
     *   jetLowC              (int)                = 0,
     *   jetHighA             (int)                = 0,
     *   jetHighB             (int)                = 0,
     *   jetHighC             (int)                = 0,
     *   gammaLowA            (int)                = 0,
     *   gammaLowB            (int)                = 0,
     *   gammaLowC            (int)                = 0,
     *   gammaHighA           (int)                = 0,
     *   gammaHighB           (int)                = 0,
     *   gammaHighC           (int)                = 0,
     *   doQA                 (bool)               = kFALSE
     */
    AliEmcalTriggerMaker *emcalTriggers = AddTaskEmcalTriggerMaker("EmcalTriggers", "EmcalTriggerSetup", 0, 0, "AliEmcalTriggerMaker", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, useOldBitConfig, doTriggerQA);
    emcalTriggers->SelectCollisionCandidates(pSel);
  }


  // ========== FILTERING ==========
//  TString dType("ESD");
//  if (!evhand->InheritsFrom("AliESDInputHandler")) 
//    dType = "AOD";

  // Set trackcuts according to period. 
  TString period(periodstr);

  TString inputClus   = TString::Format(usedClusters);
  TString inputTracks = "AODFilterTracks";
  AliEmcalAodTrackFilterTask *aodfilter =0x0;

  const Double_t edist = 440;
  if (isESD) {
    inputTracks = "ESDFilterTracks";
    TString trackCuts(Form("Hybrid_%s", period.Data()));
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalEsdTrackFilter.C");
    AliEmcalEsdTrackFilterTask *esdfilter = AddTaskEmcalEsdTrackFilter(inputTracks, trackCuts);
    esdfilter->SetDoPropagation(kTRUE);
    esdfilter->SetDist(edist);
    esdfilter->SelectCollisionCandidates(pSel);
    esdfilter->SetTrackEfficiency(trackeff);
  }
  else {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalAodTrackFilter.C");
    aodfilter = AddTaskEmcalAodTrackFilter(inputTracks, "tracks", period);
    //period here is to set specific bits for hybrid tracks. In case of HF tasks, we set below the specific bits SetAODfilterBits(16,512);
    if (doAODTrackProp) {
      aodfilter->SetDist(edist);
      aodfilter->SetAttemptPropMatch(kTRUE);
    }
    aodfilter->SelectCollisionCandidates(pSel);
    aodfilter->SetTrackEfficiency(trackeff);
    aodfilter->SetAODfilterBits(16, 512);
    aodfilter->SetMC(isMC);

  }


//----------------------- Produce EmcalParticles -----------------------------------------------------
  // // Produce objects (AliEmcalParticle) for tracks and clusters 
  // // used for cluster-track matching


  // TString emctracks = Form("EmcalTracks_%s",inputTracks.Data());
  // TString emcclusters = Form("EmcalClusters_%s",inputClus.Data());
  // Printf("emctracks: %s  inputTracks: %s",emctracks.Data(),inputTracks.Data());
  // gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalParticleMaker.C");
  // AliEmcalParticleMaker *emcalParts = AddTaskEmcalParticleMaker(inputTracks,inputClus,emctracks.Data(),emcclusters.Data());
  // emcalParts->SelectCollisionCandidates(pSel);
  // emcalParts->SetNCentBins(nCentBins);
  // // //----------------------- Cluster-Track matching -----------------------------------------------------
  // gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
  // AliEmcalClusTrackMatcherTask *emcalClus =  AddTaskEmcalClusTrackMatcher(emctracks,emcclusters,0.1,modifyMatchObjs,doHistos);
  // emcalClus->SelectCollisionCandidates(pSel);
  // emcalClus->SetNCentBins(nCentBins);

  // Printf("1-- inputTracks: %s emctracks: %s emcclusters: %s",inputTracks.Data(),emctracks.Data(),emcclusters.Data());




  
  // //----------------------- Track Matching tasks -----------------------------------------------------
  // gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskMatchingChain.C");
  // AliEmcalClusTrackMatcherTask *emcalClus =  AddTaskMatchingChain(periodstr,pSel,
  // 								  clusterColName,
  // 								  trackeff,doAODTrackProp,
  // 								  0.1,modifyMatchObjs,doHistos,nCentBins);
  
  // //hard coded names of AliEmcalParticle strings to coincide with AddTaskClusTrackMatching
  TString inputTracks = "AODFilterTracks";
  TString clusterColName = (isESD) ? "caloClusters" : "CaloClusters";
  // if (dType == "ESD") inputTracks = "ESDFilterTracks";
   TString emctracks = Form("EmcalTracks_%s", inputTracks.Data());
   TString emcclusters = Form("EmcalClusters_%s", clusterColName.Data());
  // Printf("1-- inputTracks: %s, emcclusters: %s, emctracks: %s",inputTracks.Data(),emcclusters.Data(),emctracks.Data());
  
  if (makePicoTracks) {
    TString picoTracksName(pTracksName);
    //----------------------- Produce PicoTracks -----------------------------------------------------
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
    AliEmcalPicoTrackMaker *pTrackTask = AddTaskEmcalPicoTrackMaker(picoTracksName, inputTracks);
    //    pTrackTask->SetTrackEfficiency(trackeff); //now done in Esd/AodFilter
    pTrackTask->SelectCollisionCandidates(pSel);
  }

  //----------------------- Hadronic Correction -----------------------------------------------------
  // gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskHadCorr.C"); 
  // AliHadCorrTask *hCorr = AddTaskHadCorr(emctracks,emcclusters,outClusName,hadcorr,
  // 					 minPtEt,phiMatch,etaMatch,Eexcl,trackclus,doHistos);
  // hCorr->SelectCollisionCandidates(pSel);
  // hCorr->SetNCentBins(nCentBins);

  if (isEmcalTrain) {
    if (doHistos)
      RequestMemory(hCorr, 500*1024);
  }
  
  // Produce MC particles
  TString particleColName(usedMCParticles);
  if(particleColName != "") {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskMCTrackSelector.C");
    AliEmcalMCTrackSelector *mcPartTask = AddTaskMCTrackSelector(particleColName, kFALSE, kFALSE);
    mcPartTask->SelectCollisionCandidates(pSel);
  }
  
  // Return one task that represents the jet preparation on LEGO trains
  return aodfilter;
}
