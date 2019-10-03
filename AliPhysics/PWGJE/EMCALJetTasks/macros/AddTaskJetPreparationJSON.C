// $Id$

AliAnalysisTaskSE* AddTaskJetPreparation(
  const char*    periodstr          = "LHC11h",
  const char*    pTracksName        = "PicoTracks",
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
  Bool_t   doAODTrackProp           = kFALSE,
  Bool_t   modifyMatchObjs          = kTRUE,
  Bool_t   doTriggerQA              = kFALSE
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

  //----------------------- Trigger Maker -----------------------------------------------------
  if (makeTrigger) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalTriggerMakerJSON.C");
    AliEMCALConfiguration emctriggerconf("triggerMakerConf");
    emctriggerconf.AddParam("doQA", new AliJSONBool(doTriggerQA));
    AliEmcalTriggerMaker *emcalTriggers = AddTaskEmcalTriggerMakerJSON(emctriggerconf.CreateJSONString());
    emcalTriggers->SelectCollisionCandidates(pSel);
  }

  //----------------------- Track Matching tasks -----------------------------------------------------
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskMatchingChain.C");
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
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
    AliEmcalPicoTrackMaker *pTrackTask = AddTaskEmcalPicoTrackMaker(picoTracksName, inputTracks);
    //    pTrackTask->SetTrackEfficiency(trackeff); //now done in Esd/AodFilter
    pTrackTask->SelectCollisionCandidates(pSel);
  }

  //----------------------- Hadronic Correction -----------------------------------------------------
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskHadCorr.C"); 
  AliHadCorrTask *hCorr = AddTaskHadCorr(emctracks,emcclusters,outClusName,hadcorr,
					 minPtEt,phiMatch,etaMatch,Eexcl,trackclus,doHistos);
  hCorr->SelectCollisionCandidates(pSel);
  if (isEmcalTrain) {
    if (doHistos)
      RequestMemory(hCorr,500*1024);
  }

  // Produce MC particles
  if(particleColName != "") {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskMCTrackSelector.C");
    AliEmcalMCTrackSelector *mcPartTask = AddTaskMCTrackSelector(particleColName, kFALSE, kFALSE);
    mcPartTask->SelectCollisionCandidates(pSel);
  }

  // Return one task that represents the jet preparation on LEGO trains
  return hCorr;
}
