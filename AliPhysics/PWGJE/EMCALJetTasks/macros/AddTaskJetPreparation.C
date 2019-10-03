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
  Bool_t   useOldBitConfig          = kFALSE,
  Bool_t   doTriggerQA              = kFALSE,
  Int_t    nCentBins                = 4
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
#ifdef __CLING__
    std::stringstream triggermakeradd;
    triggermakeradd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWG/EMCAL/macros/AddTaskEmcalTriggerMaker.C(";
    triggermakeradd << "\"EmcalTriggers\", \"EmcalTriggerSetup\", 0, 0, \"AliEmcalTriggerMaker\", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ";
    triggermakeradd << (useOldBitConfig ? "kTRUE" : "kFALSE") << ", ";
    triggermakeradd << (doTriggerQA ? "kTRUE" : "kFALSE");
    triggermakeradd << ")";
    std::string triggermakeraddstring = triggermakeradd.str();
    AliEmcalTriggerMaker *emcalTriggers = (AliEmcalTriggerMaker *) gROOT->ProcessLine(triggermakeraddstring.c_str());
#else
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalTriggerMaker.C");
    AliEmcalTriggerMaker *emcalTriggers = AddTaskEmcalTriggerMaker("EmcalTriggers", "EmcalTriggerSetup", 0, 0, "AliEmcalTriggerMaker", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, useOldBitConfig, doTriggerQA);
#endif
    emcalTriggers->SelectCollisionCandidates(pSel);
  }

  //----------------------- Track Matching tasks -----------------------------------------------------
#ifdef __CLING__
  std::stringstream matchingchainadd;
  matchingchainadd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWG/EMCAL/macros/AddTaskMatchingChain.C("
      << "\"" << periodstr << "\", " << pSel << ", \"" << clusterColName << "\", " << trackeff << ", "
      << (doAODTrackProp ? "kTRUE" : "kFALSE") << ", 0.1, " << (modifyMatchObjs ? "kTRUE" : "kFALSE") << ", "
      << (doHistos ? "kTRUE" : "kFALSE") << ", " << nCentBins << ")";
  std::string matchingchainaddstring = matchingchainadd.str();
  AliEmcalClusTrackMatcherTask *emcalClus = (AliEmcalClusTrackMatcherTask *)gROOT->ProcessLine(matchingchainaddstring.c_str());
#else
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskMatchingChain.C");
  AliEmcalClusTrackMatcherTask *emcalClus =  AddTaskMatchingChain(periodstr,pSel,
								  clusterColName,
								  trackeff,doAODTrackProp,
								  0.1,modifyMatchObjs,doHistos,nCentBins);
#endif
  TString inputTracks = "AODFilterTracks";
  if (dType == "ESD") inputTracks = "ESDFilterTracks";
  
  // Produce MC particles
  if(!particleColName.IsNull()) {
#ifdef __CLING__
    std::stringstream mcparttaskadd;
    mcparttaskadd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWG/EMCAL/macros/AddTaskMCTrackSelector.C("
        << "\"" << particleColName << "\", kFALSE, kFALSE)";
    std::string mcparttaskaddstring = mcparttaskadd.str();
    AliEmcalMCTrackSelector *mcPartTask = (AliEmcalMCTrackSelector *)gROOT->ProcessLine(mcparttaskaddstring.c_str());
#else
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskMCTrackSelector.C");
    AliEmcalMCTrackSelector *mcPartTask = AddTaskMCTrackSelector(particleColName, kFALSE, kFALSE);
#endif
    mcPartTask->SelectCollisionCandidates(pSel);
  }

  
  if(makePicoTracks) {
    //----------------------- Produce PicoTracks -----------------------------------------------------
#ifdef __CLING__
    std::stringstream picotrackmakeradd;
    picotrackmakeradd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C("
        << "\"" << picoTracksName << "\", \"" << inputTracks << "\")";
    std::string picotrackmakeraddstring = picotrackmakeradd.str();
    AliEmcalPicoTrackMaker *pTrackTask = (AliEmcalPicoTrackMaker *)gROOT->ProcessLine(picotrackmakeraddstring.c_str());
#else
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
    AliEmcalPicoTrackMaker *pTrackTask = AddTaskEmcalPicoTrackMaker(picoTracksName, inputTracks);
#endif
    //    pTrackTask->SetTrackEfficiency(trackeff); //now done in Esd/AodFilter
    pTrackTask->SelectCollisionCandidates(pSel);
    if(!particleColName.IsNull()) pTrackTask->SetCopyMCFlag(kTRUE, usedMCParticles);
  }

  //----------------------- Hadronic Correction -----------------------------------------------------
#ifdef __CLING__
  std::stringstream hadcorradd;
  hadcorradd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWG/EMCAL/macros/AddTaskHadCorr.C("
      << "\"" << inputTracks << "\", \"" << clusterColName << "\", \"" << outClusName << "\", " << hadcorr << ", " << minPtEt << ", "
      << phiMatch << ", " << etaMatch << ", " << Eexcl << ", " << (trackclus ? "kTRUE" : "kFALSE") << ", " << (doHistos ? "kTRUE" : "kFALSE") << ")";
  std::string hadcorraddstring = hadcorradd.str();
  AliHadCorrTask *hCorr = (AliHadCorrTask *)gROOT->ProcessLine(hadcorraddstring.c_str());
#else
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskHadCorr.C");
  AliHadCorrTask *hCorr = AddTaskHadCorr(inputTracks,clusterColName,outClusName,hadcorr,
                                         minPtEt,phiMatch,etaMatch,Eexcl,trackclus,doHistos);
#endif
  hCorr->SelectCollisionCandidates(pSel);
  hCorr->SetNCentBins(nCentBins);
#ifndef __CLING__
  if (isEmcalTrain) {
    if (doHistos)
      RequestMemory(hCorr,500*1024);
  }
#endif
  
  // Return one task that represents the jet preparation on LEGO trains
  return hCorr;
}
