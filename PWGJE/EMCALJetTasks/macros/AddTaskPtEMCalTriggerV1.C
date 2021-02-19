/**
 * \file AddTaskPtEMCalTriggerV1
 * \brief Add macro for the analysis task of high-pt tracks in triggered events
 *
 * This macro configures the analysis of high-p_{t} tracks in triggered events and adds it to
 * the analysis manager. The configuration is steered by the main function of the macro,
 * AddTaskPtEMCalTriggerV1. All other functions are helper functions performing several configuration,
 * i.e. of analysis components, and not intended for general usage.
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Dec 12, 2014
 */
#if !defined (__CINT__) || defined (__CLING__) || defined (__MAKECINT__) || defined (__ROOTCLING__)
// CLING requires header files for the JIT-compiler
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliEMCalPtTaskVTrackSelection.h"
#include "AliEMCalPtTaskTrackSelectionAOD.h"
#include "AliEMCalPtTaskTrackSelectionESD.h"
#include "AliEMCalTriggerAnaTriggerClass.h"
#include "AliEMCalTriggerAnaTriggerDecision.h"
#include "AliEMCalTriggerClusterAnalysisComponent.h"
#include "AliEMCalTriggerPatchAnalysisComponent.h"
#include "AliEMCalTriggerRecTrackAnalysisComponent.h"
#include "AliEMCalTriggerRecJetAnalysisComponent.h"
#include "AliEMCalTriggerEventCounterAnalysisComponent.h"
#include "AliEMCalTriggerMCParticleAnalysisComponent.h"
#include "AliEMCalTriggerMCJetAnalysisComponent.h"
#include "AliEMCalTriggerTaskGroup.h"
#include "AliEMCalTriggerExtraCuts.h"
#include "AliEMCalTriggerKineCuts.h"
#include "AliAnalysisTaskPtEMCalTriggerV1.h"
#include "AliESDtrackCuts.h"
#include "AliClusterContainer.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliVEvent.h"
#include <TList.h>
#include <TString.h>
#include <cstring>
#endif

void AddClusterComponent(PWGJE::EMCALJetTasks::AliEMCalTriggerTaskGroup *group);
void AddTrackComponent(PWGJE::EMCALJetTasks::AliEMCalTriggerTaskGroup *group, AliEmcalTrackSelection *trackcuts, bool isMC, bool isSwapEta);
void AddMCParticleComponent(PWGJE::EMCALJetTasks::AliEMCalTriggerTaskGroup *group);
void AddEventCounterComponent(PWGJE::EMCALJetTasks::AliEMCalTriggerTaskGroup *group);
void AddMCJetComponent(PWGJE::EMCALJetTasks::AliEMCalTriggerTaskGroup *group, double minJetPt);
void AddRecJetComponent(PWGJE::EMCALJetTasks::AliEMCalTriggerTaskGroup *group, AliEmcalTrackSelection *trackcuts, double minJetPt, bool isMC, bool isSwapEta);
void CreateJetPtBinning(PWGJE::EMCALJetTasks::AliAnalysisTaskPtEMCalTriggerV1 *task);
void CreateTriggerClassespPb2013(PWGJE::EMCALJetTasks::AliAnalysisTaskPtEMCalTriggerV1 *task, bool isMC);
void CreateTriggerClassespp2012(PWGJE::EMCALJetTasks::AliAnalysisTaskPtEMCalTriggerV1 *task, bool isMC);
AliEmcalTrackSelection *CreateDefaultTrackCuts(bool isAOD);
AliEmcalTrackSelection *CreateHybridTrackCuts(bool isAOD);
AliEmcalTrackSelection *TrackCutsFactory(const char *trackCutsName, bool isAOD);

/**
 * \brief Configuring the analysis of high-p_{t} tracks in triggered events and adds it to the analysis train
 *
 * Main function of the macro: Configures the task (enabling components, linking connections, defining output containers)
 * and adds it to the analysis manager. The macro auto-configures the task for ESD or AOD analysis. The user has to define
 * only whether the analysis runs on Monte-Carlo, and whether the Monte-Carlo sample is done in p_{t}-hard bins. Sample-
 * specific settings are steered via the string period.
 *
 * The function enables several analysis components. These are:
 *  - particles: Analysis of MC true particles
 *  - triggers: Analysis of trigger patches
 *  - clusters: Analysis of EMCal clusters
 *  - tracks: Analysis of reconstructed tracks
 *  - mcjets: Analysis of jets build from MC-true particles
 *  - recjets: Analysis of jet build from reconstructed tracks
 *  The string of components (see under parameters) uses : as separator. All analysis components are linked to task groups which
 *  share a common event selection. The task contains two groups, one without event selection and a second with a standard event
 *  selection (mainly vertex selection). The trigger patch analysis component, if enabled, is connected to the group without any event
 *  selection, while all other analyses components are connected to the group with the default event selection. By default all
 *  components are switched on.
 *
 * \param isMC True if MC-truth is available
 * \param usePythiaHard True in case the analysis is done on a production in pt-hard bins
 * \param period Name of the period
 * \param ntrackContainer Name of the track container
 * \param nclusterContainer Name of the cluster container
 * \param njetcontainerData Name of the jet container on reconstructed event
 * \param njetcontainerMC Name of the jet container on true event
 * \param ntriggerContainer Name of the container having trigger patches
 * \param jetradius Jet resolution parameter used for both jetfinders
 * \param ntrackcuts Name of the track cuts used in the selection of reconstructed tracks
 * \param components Components enabled in the analysis (see above)
 * \param triggersetup Pre-defined trigger setup for several periods
 * \param useOfflinePatches True if offline patches are used (only relevant in case patches are used for the trigger decision)
 * \return The fully configured task
 */
AliAnalysisTask* AddTaskPtEMCalTriggerV1(
    bool isMC,
    bool usePythiaHard,
    const char *period ="LHC13d",
    const char *ntrackContainer = "",
    const char *nclusterContainer = "",
    const char *njetcontainerData = "",
    const char *njetcontainerMC = "",
    const char *ntriggerContainer = "",
    double jetradius = 0.5,
    const char *ntrackcuts = "standard",
    const char *components = "particles:clusters:tracks:mcjets:recjets:triggers:patchevent",
    const char * triggersetup = "pPb2013",
    bool useOfflinePatches = kFALSE,
    const char *options = ""
)
{
  //AliLog::SetClassDebugLevel("PWGJE::EMCALJetTasks::AliAnalysisTaskPtEMCalTrigger", 2);
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    ::Error("AddTaskPtEMCalTrigger", "No analysis manager to connect to.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPtEMCalTrigger", "This task requires an input event handler");
    return NULL;
  }
  bool isAOD = (mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class());

  // Decode components
  bool doClusters(false), doMCParticles(false), doTracks(false), doMCJets(false), doRecJets(false), doTriggers(false), doTriggersEvents(false);
  TObjArray *compsplit = TString(components).Tokenize(":");
  TIter tokenIter(compsplit);
  TObjString *compstring(NULL);
  while((compstring = (TObjString *)tokenIter())){
	  if(!compstring->String().CompareTo("clusters")) doClusters = true;
	  if(!compstring->String().CompareTo("particles")) doMCParticles = true;
	  if(!compstring->String().CompareTo("tracks")) doTracks = true;
	  if(!compstring->String().CompareTo("mcjets")) doMCJets = true;
	  if(!compstring->String().CompareTo("recjets")) doRecJets = true;
	  if(!compstring->String().CompareTo("triggers")) doTriggers = true;
	  if(!compstring->String().CompareTo("patchevent")) doTriggersEvents = true;
  }

  std::cout << "Track task configuration:" << std::endl;
  std::cout << "==================================" << std::endl;
  std::cout << "Monte-Carlo particles:      " << (doMCParticles ? "ON" : "OFF") << std::endl;
  std::cout << "Monte-Carlo jets:           " << (doMCJets ? "ON" : "OFF") << std::endl;
  std::cout << "Trigger Patch QA:           " << (doTriggers ? "ON" : "OFF") << std::endl;
  std::cout << "Tracks:                     " << (doTracks ? "ON" : "OFF") << std::endl;
  std::cout << "EMCAL clusters:             " << (doClusters ? "ON" : "OFF") << std::endl;
  std::cout << "Reconstructed jets:         " << (doRecJets ? "ON" : "OFF") << std::endl;

  bool isSwapEta = TString(period).CompareTo("LHC13f") ? kFALSE : kTRUE;
  PWGJE::EMCALJetTasks::AliAnalysisTaskPtEMCalTriggerV1 *pttriggertask = new PWGJE::EMCALJetTasks::AliAnalysisTaskPtEMCalTriggerV1(Form("ptemcaltriggertask%s", ntrackcuts));
  //pttriggertask->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kEMC7);                          // Select both INT7 or EMC7 triggered events
  pttriggertask->SelectCollisionCandidates(AliVEvent::kAny);

  PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerDecisionConfig *trgconf = new PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerDecisionConfig;
  if(isMC && !useOfflinePatches) trgconf->SetSwapThresholds();
  //printf("Using offline patches: %s\n", useOfflinePatches ? "yes" :"no");
  trgconf->SetUseOfflinePatches(useOfflinePatches);
  if(!useOfflinePatches){
    // Cut further in amplitude calibrated amplituded on the trigger patches
    trgconf->SetPatchEnergyType(PWGJE::EMCALJetTasks::kAmplitudeOffline);
    trgconf->SetEnergyThreshold(PWGJE::EMCALJetTasks::kTAEMCGHigh, 140);
    trgconf->SetEnergyThreshold(PWGJE::EMCALJetTasks::kTAEMCGLow, 89);
    trgconf->SetEnergyThreshold(PWGJE::EMCALJetTasks::kTAEMCJHigh, 260);
    trgconf->SetEnergyThreshold(PWGJE::EMCALJetTasks::kTAEMCJLow, 127);
  }
  pttriggertask->SetTriggerDecisionConfig(trgconf);

  // Create trigger setup
  if(!TString(triggersetup).CompareTo("pPb2013")){
   CreateTriggerClassespPb2013(pttriggertask, isMC);
  } else if(!TString(triggersetup).CompareTo("pp2012")){
   CreateTriggerClassespp2012(pttriggertask, isMC);
  } else {
    std::cout << "No trigger setup defined for the given request " << triggersetup << std::endl;
  }

  CreateJetPtBinning(pttriggertask);
//  pttriggertask->SetTriggerDebug(kTRUE);

  mgr->AddTask(pttriggertask);
  if(usePythiaHard){
    pttriggertask->SetIsPythia(kTRUE);
  }

  if(strlen(ntriggerContainer)){
    pttriggertask->SetCaloTriggerPatchInfoName(ntriggerContainer);
  }

  // Add components
  if(doTriggers){
    PWGJE::EMCALJetTasks::AliEMCalTriggerTaskGroup *noselect = new PWGJE::EMCALJetTasks::AliEMCalTriggerTaskGroup("noselect");
    PWGJE::EMCALJetTasks::AliEMCalTriggerPatchAnalysisComponent * triggercomp = new PWGJE::EMCALJetTasks::AliEMCalTriggerPatchAnalysisComponent("patchanalysis");
    if(isMC) triggercomp->SetSwapOnlineThresholds(true);
    noselect->AddAnalysisComponent(triggercomp);
    pttriggertask->AddAnalysisGroup(noselect);
  }

  double jetpt[4] = {40., 60., 80., 100.};
  PWGJE::EMCALJetTasks::AliEMCalTriggerTaskGroup *defaultselect = new PWGJE::EMCALJetTasks::AliEMCalTriggerTaskGroup("defaultselect");
  PWGJE::EMCALJetTasks::AliEMCalTriggerEventSelection *eventselect = PWGJE::EMCALJetTasks::AliEMCalTriggerEventSelection();
  TString optstring(options);
  if(optstring.Contains("oldpileup")) eventselect->SetOldPileupSelection();
  if(optstring.Contains("oldvertex")) eventselect->SetOldVertexSelection();
  defaultselect->SetEventSelection(eventselect);
  PWGJE::EMCALJetTasks::AliEMCalTriggerKineCuts *kineCuts = new PWGJE::EMCALJetTasks::AliEMCalTriggerKineCuts();
  kineCuts->SetPtRange(2., 100.);
  defaultselect->SetKineCuts(kineCuts);
  AddEventCounterComponent(defaultselect);
  if(isMC){
    if(doMCParticles) AddMCParticleComponent(defaultselect);
    if(doMCJets) AddMCJetComponent(defaultselect, 20.);
    /*
    for(int ijpt = 0; ijpt < 4; ijpt++)
      AddMCJetComponent(defaultselect, jetpt[ijpt]);
    */
  }
  if(doClusters) 	AddClusterComponent(defaultselect);
  if(doTracks) 		AddTrackComponent(defaultselect, TrackCutsFactory(ntrackcuts, isAOD), isMC, isSwapEta);
  if(doRecJets) 	AddRecJetComponent(defaultselect, TrackCutsFactory(ntrackcuts, isAOD), 20., isMC, isSwapEta);
  /*
   * for(int ijpt = 0; ijpt < 4; ijpt++)
       AddRecJetComponent(defaultselect, TrackCutsFactory(ntrackcuts), jetpt[ijpt], isMC, isSwapEta);
   */
  if(doTriggersEvents){
    PWGJE::EMCALJetTasks::AliEMCalTriggerPatchAnalysisComponent * triggereventcomp = new PWGJE::EMCALJetTasks::AliEMCalTriggerPatchAnalysisComponent("patchineventanalysis", kTRUE);
    if(isMC) triggereventcomp->SetSwapOnlineThresholds(true);
    defaultselect->AddAnalysisComponent(triggereventcomp);
  }

  pttriggertask->AddAnalysisGroup(defaultselect);

  // Add containers
  Bool_t isAOD = mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class();
  AliParticleContainer *trackContainer = pttriggertask->AddParticleContainer(ntrackContainer);
  //trackContainer->SetClassName("AliVTrack");
  AliClusterContainer *clusterContainer = pttriggertask->AddClusterContainer(nclusterContainer);
  AliParticleContainer *mcpartcont = isMC ? pttriggertask->AddParticleContainer("MCParticlesSelected") : NULL;


  // Handle Jet Containers
  if(strlen(njetcontainerData)){
    AliJetContainer *jetcontainerData = pttriggertask->AddJetContainer(njetcontainerData, "EMCAL", jetradius);
    pttriggertask->SetDataJetContainerName("PtTriggerTaskJetsData");
    jetcontainerData->ConnectParticleContainer(trackContainer);
    jetcontainerData->SetName("PtTriggerTaskJetsData");
    jetcontainerData->SetJetPtCut(20.);
    printf("jet container added for Data\n");
  }
  if(isMC && strlen(njetcontainerMC)){
    AliJetContainer *jetcontainerMC = pttriggertask->AddJetContainer(njetcontainerMC, "EMCAL", jetradius);
    pttriggertask->SetMCJetContainerName("PtTriggerTaskJetsMC");
    jetcontainerMC->ConnectParticleContainer(mcpartcont);
    jetcontainerMC->SetName("PtTriggerTaskJetsMC");
    jetcontainerMC->SetJetPtCut(20.);
    printf("Jet container added for MC");
  }

  TString containerName = mgr->GetCommonFileName();
  containerName += ":PtEMCalTriggerTask" + TString(ntrackcuts);
  printf("container name: %s\n", containerName.Data());

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("TriggerTracksResults%s", ntrackcuts), TList::Class(),    AliAnalysisManager::kOutputContainer, containerName.Data());

  //Connect input/output
  mgr->ConnectInput(pttriggertask, 0, cinput);
  mgr->ConnectOutput(pttriggertask, 1, coutput);

  return pttriggertask;
}

/**
 * \brief Configuration of analysis component on EMCAL clusters
 *
 * Configures the analysis component on EMCAL clusters (calibrated and uncalibrated) and adds it to the parent analysis group. Function
 * is a helper function called by AddTaskPtEMCalTriggerV1 and not intended for general usage.
 *
 * \param group Parent group the component is linked to
 */
void AddClusterComponent(PWGJE::EMCALJetTasks::AliEMCalTriggerTaskGroup *group){
  PWGJE::EMCALJetTasks::AliEMCalTriggerClusterAnalysisComponent *clusteranalysis = new PWGJE::EMCALJetTasks::AliEMCalTriggerClusterAnalysisComponent("clusterAnalysis");
  clusteranalysis->SetEnergyRange(2., 100.);
  group->AddAnalysisComponent(clusteranalysis);
}

/**
 * \brief Configuration of analysis component on reconstructed tracks
 *
 * Configures the analysis component on reconstructed tracks and adds it to its parent analysis group. Function
 * is a helper function called by AddTaskPtEMCalTriggerV1 and not intended for general usage.
 *
 * \param group Parent analysis group
 * \param trackcuts Cuts used for the track selection
 * \param isMC True if MC information is available.
 * \param isSwapEta True if the \f$ \eta \f$ sign is swapped
 */
void AddTrackComponent(PWGJE::EMCALJetTasks::AliEMCalTriggerTaskGroup *group, AliEmcalTrackSelection * trackcuts, bool isMC, bool isSwapEta){
  PWGJE::EMCALJetTasks::AliEMCalTriggerRecTrackAnalysisComponent *trackanalysis = new PWGJE::EMCALJetTasks::AliEMCalTriggerRecTrackAnalysisComponent("trackAnalysisStandard");
  group->AddAnalysisComponent(trackanalysis);
  // Create charged hadrons pPb standard track cuts
  trackanalysis->SetTrackSelection(trackcuts);

  if(isMC) trackanalysis->SetRequestMCtrueTracks();
  if(isSwapEta) trackanalysis->SetSwapEta();
}

/**
 * \brief Configuration of event counter component.
 *
 * Configures event counter component and adds it to the parent analysis group. Function is a helper function
 * called by AddTaskPtEMCalTriggerV1 and not intended for general usage.
 *
 * \param group Parent group of analysis components
 */
void AddEventCounterComponent(PWGJE::EMCALJetTasks::AliEMCalTriggerTaskGroup *group){
  PWGJE::EMCALJetTasks::AliEMCalTriggerEventCounterAnalysisComponent * evcount = new PWGJE::EMCALJetTasks::AliEMCalTriggerEventCounterAnalysisComponent("eventCounter");
  group->AddAnalysisComponent(evcount);
}

/**
 * \brief Configuration of analysis component on MC true events and particles
 *
 * Configures the analysis component on MC true events and particles and adds it to the parent analysis group. Function is a
 * helper function called by AddTaskPtEMCalTriggerV1 and not intended for general usage.
 *
 * \param group The parent analysis group
 */
void AddMCParticleComponent(PWGJE::EMCALJetTasks::AliEMCalTriggerTaskGroup *group){
  group->AddAnalysisComponent(new PWGJE::EMCALJetTasks::AliEMCalTriggerMCParticleAnalysisComponent("partana"));
}

/**
 * \brief Configuration of analysis component on MC true jets
 *
 * Configures the analysis component on MC true jets (jets from MC true particles in MC true events) and adds it to the analysis group.
 * Function is a helper function called by AddTaskPtEMCalTriggerV1 and not intended for general usage.
 *
 * \param group The parent analysis group
 * \param minJetPt Lower jet-p_{t} cut of the analysis
 */
void AddMCJetComponent(PWGJE::EMCALJetTasks::AliEMCalTriggerTaskGroup *group, double minJetPt){
  PWGJE::EMCALJetTasks::AliEMCalTriggerMCJetAnalysisComponent *jetana = new PWGJE::EMCALJetTasks::AliEMCalTriggerMCJetAnalysisComponent(Form("MCJetAna%f", minJetPt));
  jetana->SetMinimumJetPt(minJetPt);
  group->AddAnalysisComponent(jetana);
}

/**
 * \brief Configuration of analysis on reconstructed jets
 *
 * Configures the analysis on reconstructed jets and adds it to the parent analysis group. Function is a
 * helper function called by AddTaskPtEMCalTriggerV1 and not intended for general usage.
 *
 * \param group The parent analysis group
 * \param trackcuts Track selection applied to particles in reconstructed jets
 * \param minJetPt Lower jet-p_{t} cut of the analysis
 * \param isMC True if MC information is available.
 * \param isSwapEta True if the \f$ \eta \f$ sign is swapped
 */
void AddRecJetComponent(PWGJE::EMCALJetTasks::AliEMCalTriggerTaskGroup *group, AliEmcalTrackSelection *trackcuts, double minJetPt, bool isMC, bool isSwapEta){
  PWGJE::EMCALJetTasks::AliEMCalTriggerRecJetAnalysisComponent *jetana = new PWGJE::EMCALJetTasks::AliEMCalTriggerRecJetAnalysisComponent(Form("RecJetAna%f", minJetPt));
  jetana->SetMinimumJetPt(minJetPt);
  jetana->SetSingleTrackCuts(trackcuts);
  //jetana->SetComponentDebugLevel(2);
  group->AddAnalysisComponent(jetana);
}

/**
 * \brief Create binning for jet p_{t} axis
 *
 * Create binning for jet p_{t} axis and adds it to the binning factory of the parent task. Function is a
 * helper function called by AddTaskPtEMCalTriggerV1 and not intended for general usage.
 *
 * \param task Task the binning is associated to
 */
void CreateJetPtBinning(PWGJE::EMCALJetTasks::AliAnalysisTaskPtEMCalTriggerV1 *task){
  // Linear binning in steps of 10 GeV/c up to 200 GeV/c
  TArrayD binlimits(21);
  for(int i = 0; i < 21; i++) binlimits[i] = 10.*i;
  task->SetBinning("jetpt", binlimits);
}

/**
 * Define trigger class handling for p-Pb 2013 and add the classes to the analysis task. Trigger classes available:
 *  - 1 Min bias trigger
 *  - EMCAL jet trigger, high threshold
 *  - EMCAL jet trigger, low threshold
 *  - EMCAL single shower trigger, high threshold
 *  - EMCAL single shower trigger, low threshold
 *  In case of data the EMCAL trigger decision is obtained from the trigger string while in case of MC it is obtained
 *  from the trigger patches.
 *
 * \param task Target task
 * \param isMC Definition whether task runs on data or on MC
 */
void CreateTriggerClassespPb2013(PWGJE::EMCALJetTasks::AliAnalysisTaskPtEMCalTriggerV1 *task, bool isMC){
  // Min Bias trigger
  PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass *minbiastrigger = new PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass("MinBias", "min. bias events");
  minbiastrigger->SetMinBiasTrigger(kTRUE);
  minbiastrigger->AddTriggerBit(AliVEvent::kINT7);
  task->AddTriggerClass(minbiastrigger);

  // EMCAL trigger classes
  PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass *jethigh = new PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass("EMCJHigh", "jet-triggered events (high threshold)");
  PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass *jetlow = new PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass("EMCJLow", "jet-triggered events (low threshold)");
  PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass *gammahigh = new PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass("EMCGHigh", "gamma-triggered events (high threshold)");
  PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass *gammalow = new PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass("EMCGLow", "gamma-triggered events (low threshold)");

  if(isMC){
    // Get triggers from patches
    jethigh->AddTriggerStringPattern("EJ1", kTRUE);
    jetlow->AddTriggerStringPattern("EJ2", kTRUE);
    gammahigh->AddTriggerStringPattern("EG1", kTRUE);
    gammalow->AddTriggerStringPattern("EG2", kTRUE);
  } else {
    // Get triggers from trigger string
    jethigh->AddTriggerPatchType(PWGJE::EMCALJetTasks::kTAEMCJHigh);
    jethigh->AddTriggerPatchType(PWGJE::EMCALJetTasks::kTAEMCJLow);
    jethigh->AddTriggerPatchType(PWGJE::EMCALJetTasks::kTAEMCGHigh);
    jethigh->AddTriggerPatchType(PWGJE::EMCALJetTasks::kTAEMCGLow);
  }
  task->AddTriggerClass(jethigh);
  task->AddTriggerClass(jetlow);
  task->AddTriggerClass(gammahigh);
  task->AddTriggerClass(gammalow);

  /* No handling for mixed classes for the moment, these would be
   *
   * "EMCHighBoth", "jet and gamma triggered events (high threshold)"
   * "EMCHighGammaOnly", "exclusively gamma-triggered events (high threshold)"
   * "EMCHighJetOnly", "exclusively jet-triggered events (high threshold)"
   * "EMCLowBoth", "jet and gamma triggered events (low threshold)"
   * "EMCLowGammaOnly", "exclusively gamma-triggered events (low threshold)"
   * "EMCLowJetOnly", "exclusively-triggered events (low threshold)"
   */
}

/**
 * Define trigger class handling for pp 8 TeV 2012 and add the classes to the analysis task. Trigger classes available:
 *  - 2 Min bias triggers (INT7 and INT8, handled separately)
 *  - EMCAL jet trigger, level 0 based on INT7 or INT8 (handled separately)
 *  - EMCAL single shower trigger, level 0 based on INT7 or INT8 (handled separately)
 *  In case of data the EMCAL trigger decision is obtained from the trigger string while in case of MC it is obtained
 *  from the trigger patches and the min. bias bit for INT7 or INT8.
 *
 * \param task Target task
 * \param isMC Definition whether task runs on data or on MC
 */
void CreateTriggerClassespp2012(PWGJE::EMCALJetTasks::AliAnalysisTaskPtEMCalTriggerV1 *task, bool isMC){
  // Min Bias triggers (INT7 and INT8)
  PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass * mbINT7 = new PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass("MinBiasINT7", "min. bias events (INT7)");
  mbINT7->SetMinBiasTrigger(kTRUE);
  mbINT7->AddTriggerBit(AliVEvent::kINT7);
  task->AddTriggerClass(mbINT7);
  PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass * mbINT8 = new PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass("MinBiasINT8", "min. bias events (INT8)");
  mbINT8->SetMinBiasTrigger(kTRUE);
  mbINT8->AddTriggerBit(AliVEvent::kINT8);
  task->AddTriggerClass(mbINT8);

  PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass * jetINT7 = new PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass("EMCJetINT7", "EMCAL jet triggered events (INT7)");
  PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass * jetINT8 = new PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass("EMCJetINT8", "EMCAL jet triggered events (INT8)");
  PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass * gammaINT7 = new PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass("EMCGammaINT7", "EMCAL gamma triggered events (INT7)");
  PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass * gammaINT8 = new PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerClass("EMCGammaINT8", "EMCAL gamma triggered events (INT8)");
  if(!isMC){
    // EMCAL triggers for single shower, separated for INT7 and INT8, detected from trigger string
    jetINT7->AddTriggerStringPattern("EJE", kTRUE);
    jetINT7->AddTriggerStringPattern("CEMC7", kTRUE);
    jetINT8->AddTriggerStringPattern("EJE", kTRUE);
    jetINT8->AddTriggerStringPattern("CEMC8", kTRUE);
    gammaINT7->AddTriggerStringPattern("EGA", kTRUE);
    gammaINT7->AddTriggerStringPattern("CEMC7", kTRUE);
    gammaINT8->AddTriggerStringPattern("EGA", kTRUE);
    gammaINT8->AddTriggerStringPattern("CEMC8", kTRUE);
  } else {
    jetINT7->AddTriggerPatchType(PWGJE::EMCALJetTasks::kTAEMCJHigh);
    jetINT7->AddTriggerBit(AliVEvent::kINT7);
    jetINT8->AddTriggerPatchType(PWGJE::EMCALJetTasks::kTAEMCJHigh);
    jetINT8->AddTriggerBit(AliVEvent::kINT8);
    gammaINT7->AddTriggerPatchType(PWGJE::EMCALJetTasks::kTAEMCGHigh);
    gammaINT7->AddTriggerBit(AliVEvent::kINT7);
    gammaINT8->AddTriggerPatchType(PWGJE::EMCALJetTasks::kTAEMCGHigh);
    gammaINT8->AddTriggerBit(AliVEvent::kINT8);
  }

  task->AddTriggerClass(jetINT7);
  task->AddTriggerClass(jetINT8);
  task->AddTriggerClass(gammaINT7);
  task->AddTriggerClass(gammaINT8);
}

/**
 * \brief Creation of default (\f$ R_{pPb} \f$) track selection cuts
 *
 * Creates the default (\f$ R_{pPb} \f$) track selection for ESD or AOD analysis. In order to be independent of the
 * analysis type, the track selection is hidden from the analysis flow inside a wrapper for track selections, which
 * is implemented for ESD and AOD analyses. Function is a  helper function called by AddTaskPtEMCalTriggerV1 and not
 * intended for general usage.
 *
 * \param isAOD True in case the analysis is performed on AOD tracks
 * \return A virtual track selection using the hybrid track selection cuts
 */
AliEmcalTrackSelection *CreateDefaultTrackCuts(bool isAOD){
  AliEmcalTrackSelection * trackSelection(NULL);
  if(isAOD){
	  AliEmcalTrackSelectionAOD *aodsel = new AliEmcalTrackSelectionAOD;
	  aodsel->AddFilterBit(AliAODTrack::kTrkGlobal);
	  PWGJE::EMCALJetTasks::AliEMCalTriggerExtraCuts *extraCuts = new PWGJE::EMCALJetTasks::AliEMCalTriggerExtraCuts();
	  extraCuts->SetMinTPCCrossedRows(120);
	  //extraCuts->SetMinTPCTrackLengthCut();
	  aodsel->AddTrackCuts(extraCuts);
	  trackSelection = aodsel;
  } else {
	  AliESDtrackCuts *standardTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
	  standardTrackCuts->SetName("Standard Track cuts");
	  standardTrackCuts->SetMinNCrossedRowsTPC(120);
	  standardTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
	  trackSelection = new AliEmcalTrackSelectionESD(standardTrackCuts);
	  //PWGJE::EMCALJetTasks::AliEMCalTriggerExtraCuts *extraCuts = new PWGJE::EMCALJetTasks::AliEMCalTriggerExtraCuts();
	  //extraCuts->SetMinTPCTrackLengthCut();
	  //trackSelection->AddTrackCuts(extraCuts);
  }
  return trackSelection;
}

/**
 * \brief Creation of default (\f$ R_{pPb} \f$) track selection cuts
 *
 * Creates the hybrid track selection for ESD or AOD analysis. In order to be independent of the analysis type,
 * the track selection is hidden from the analysis flow inside a wrapper for track selections, which is implemented
 * for ESD and AOD analyses. Function is a  helper function called by AddTaskPtEMCalTriggerV1 and not intended for
 * general usage.
 *
 * \param isAOD True in case the analysis is performed on AOD tracks
 * \return A virtual track selection using the hybrid track selection cuts
 */
AliEmcalTrackSelection *CreateHybridTrackCuts(bool isAOD){
  AliEmcalTrackSelection * trackSelection(NULL);
  if(isAOD){
	  // Purely use filter bits
	  AliEmcalTrackSelectionAOD *aodsel = new AliEmcalTrackSelectionAOD;
	  aodsel->AddFilterBit(256);
	  aodsel->AddFilterBit(512);
	  trackSelection = aodsel;
  } else {
	  AliESDtrackCuts* hybridTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
	  hybridTrackCuts->SetName("Global Hybrid tracks, loose DCA");
	  hybridTrackCuts->SetMaxDCAToVertexXY(2.4);
	  hybridTrackCuts->SetMaxDCAToVertexZ(3.2);
	  hybridTrackCuts->SetDCAToVertex2D(kTRUE);
	  hybridTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
	  hybridTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
	  trackSelection = new AliEmcalTrackSelectionESD(hybridTrackCuts);
  }
  return trackSelection;
}

/**
 * \brief Steering function of track cuts creation
 *
 * This function steers the creation of the track selection objects and delegates it to other function. Currently
 * the following track selections are implemented:
 *  - standard: The standard (\f$ R_{pPb} \f$) track selection
 *  - hybrid: The hybrid track selection, used in jet analyses
 * Function is a  helper function called by AddTaskPtEMCalTriggerV1 and not intended for general usage.
 *
 * \param trackCutsName Name of the track cuts
 * \param isAOD True in case the analysis is performed on AOD tracks
 * \return A virtual track selection object (NULL for invalid track cut names)
 */
AliEmcalTrackSelection *TrackCutsFactory(const char* trackCutsName, bool isAOD) {
  if(!strcmp(trackCutsName, "standard")) return CreateDefaultTrackCuts(isAOD);
  else if(!strcmp(trackCutsName, "hybrid")) return CreateHybridTrackCuts(isAOD);
  return NULL;
}
