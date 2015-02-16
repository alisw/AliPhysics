#if !defined (__CINT__) || defined (__MAKECINT__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskPtEMCalTriggerV1.h"
#include "AliESDtrackCuts.h"
#include "AliJetContainer.h"
#include <TList.h>
#include <TString.h>
#include <cstring>
#endif

void AddClusterComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, bool usePatches);
void AddTrackComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, EMCalTriggerPtAnalysis::AliEMCalPtTaskVTrackSelection *trackcuts, bool isMC, bool usePatches, bool isSwapEta);
void AddMCParticleComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group);
void AddEventCounterComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, bool usePatches);
void AddMCJetComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, double minJetPt);
void AddRecJetComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, EMCalTriggerPtAnalysis::AliEMCalPtTaskVTrackSelection *trackcuts, double minJetPt, bool isMC, bool usePatches, bool isSwapEta);
void CreateJetPtBinning(EMCalTriggerPtAnalysis::AliAnalysisTaskPtEMCalTriggerV1 *task);
EMCalTriggerPtAnalysis::AliEMCalPtTaskVTrackSelection *CreateDefaultTrackCuts(bool isAOD);
EMCalTriggerPtAnalysis::AliEMCalPtTaskVTrackSelection *CreateHybridTrackCuts(bool isAOD);
EMCalTriggerPtAnalysis::AliEMCalPtTaskVTrackSelection *TrackCutsFactory(const char *trackCutsName, bool isAOD);

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
	const char *components = "particles:clusters:tracks:mcjets:recjets:triggers",
	bool usePatches = kFALSE,
	bool useOfflinePatches = kFALSE
)
{
  //AliLog::SetClassDebugLevel("EMCalTriggerPtAnalysis::AliAnalysisTaskPtEMCalTrigger", 2);
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
  bool doClusters(false), doMCParticles(false), doTracks(false), doMCJets(false), doRecJets(false), doTriggers(false);
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
  EMCalTriggerPtAnalysis::AliAnalysisTaskPtEMCalTriggerV1 *pttriggertask = new EMCalTriggerPtAnalysis::AliAnalysisTaskPtEMCalTriggerV1(Form("ptemcaltriggertask%s", ntrackcuts));
  //pttriggertask->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kEMC7);                          // Select both INT7 or EMC7 triggered events
  pttriggertask->SelectCollisionCandidates(AliVEvent::kAny);

  EMCalTriggerPtAnalysis::AliEMCalTriggerAnaTriggerDecisionConfig *trgconf = new EMCalTriggerPtAnalysis::AliEMCalTriggerAnaTriggerDecisionConfig;
  if(isMC) trgconf->SetSwapThresholds();
  trgconf->SetUseOfflinePatches(useOfflinePatches);
  pttriggertask->SetTriggerDecisionConfig(trgconf);

  CreateJetPtBinning(pttriggertask);
  //pttriggertask->SetTriggerDebug(kTRUE);

  mgr->AddTask(pttriggertask);
  if(usePythiaHard){
    pttriggertask->SetIsPythia(kTRUE);
  }

  if(strlen(ntriggerContainer)){
    pttriggertask->SetCaloTriggerPatchInfoName(ntriggerContainer);
  }

  // Add components
  if(doTriggers){
    EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *noselect = new EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup("noselect");
    noselect->AddAnalysisComponent(new EMCalTriggerPtAnalysis::AliEMCalTriggerPatchAnalysisComponent("patchanalysis"));
    pttriggertask->AddAnalysisGroup(noselect);
  }

  double jetpt[4] = {40., 60., 80., 100.};
  EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *defaultselect = new EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup("defaultselect");
  defaultselect->SetEventSelection(new EMCalTriggerPtAnalysis::AliEMCalTriggerEventSelection());
  EMCalTriggerPtAnalysis::AliEMCalTriggerKineCuts *kineCuts = new EMCalTriggerPtAnalysis::AliEMCalTriggerKineCuts();
  kineCuts->SetPtRange(2., 100.);
  defaultselect->SetKineCuts(kineCuts);
  AddEventCounterComponent(defaultselect, usePatches);
  if(isMC){
    if(doMCParticles) AddMCParticleComponent(defaultselect);
    if(doMCJets) AddMCJetComponent(defaultselect, 20.);
    /*
    for(int ijpt = 0; ijpt < 4; ijpt++)
      AddMCJetComponent(defaultselect, jetpt[ijpt]);
    */
  }
  if(doClusters) 	AddClusterComponent(defaultselect, usePatches);
  if(doTracks) 		AddTrackComponent(defaultselect, TrackCutsFactory(ntrackcuts, isAOD), isMC, usePatches, isSwapEta);
  if(doRecJets) 	AddRecJetComponent(defaultselect, TrackCutsFactory(ntrackcuts, isAOD), 20., isMC, usePatches, isSwapEta);
  /*
   * for(int ijpt = 0; ijpt < 4; ijpt++)
       AddRecJetComponent(defaultselect, TrackCutsFactory(ntrackcuts), jetpt[ijpt], isMC, isSwapEta);
   */

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

void AddClusterComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, bool usePatches){
  EMCalTriggerPtAnalysis::AliEMCalTriggerClusterAnalysisComponent *clusteranalysis = new EMCalTriggerPtAnalysis::AliEMCalTriggerClusterAnalysisComponent("clusterAnalysis");
  clusteranalysis->SetEnergyRange(2., 100.);
  clusteranalysis->SetUsePatches(usePatches);
  group->AddAnalysisComponent(clusteranalysis);
}

void AddTrackComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, EMCalTriggerPtAnalysis::AliEMCalPtTaskVTrackSelection * trackcuts, bool isMC, bool usePatches, bool isSwapEta){
  EMCalTriggerPtAnalysis::AliEMCalTriggerRecTrackAnalysisComponent *trackanalysis = new EMCalTriggerPtAnalysis::AliEMCalTriggerRecTrackAnalysisComponent("trackAnalysisStandard");
  group->AddAnalysisComponent(trackanalysis);
  // Create charged hadrons pPb standard track cuts
  trackanalysis->SetTrackSelection(trackcuts);

  if(isMC) trackanalysis->SetRequestMCtrueTracks();
  if(usePatches) trackanalysis->SetUsePatches();
  if(isSwapEta) trackanalysis->SetSwapEta();
}

void AddEventCounterComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, Bool_t usePatches){
  EMCalTriggerPtAnalysis::AliEMCalTriggerEventCounterAnalysisComponent * evcount = new EMCalTriggerPtAnalysis::AliEMCalTriggerEventCounterAnalysisComponent("eventCounter");
  evcount->SetUsePatches(usePatches);
  group->AddAnalysisComponent(evcount);
}

void AddMCParticleComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group){
  group->AddAnalysisComponent(new EMCalTriggerPtAnalysis::AliEMCalTriggerMCParticleAnalysisComponent("partana"));
}

void AddMCJetComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, double minJetPt){
  EMCalTriggerPtAnalysis::AliEMCalTriggerMCJetAnalysisComponent *jetana = new EMCalTriggerPtAnalysis::AliEMCalTriggerMCJetAnalysisComponent(Form("MCJetAna%f", minJetPt));
  jetana->SetMinimumJetPt(minJetPt);
  jetana->SetUsePatches();
  group->AddAnalysisComponent(jetana);
}

void AddRecJetComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, EMCalTriggerPtAnalysis::AliEMCalPtTaskVTrackSelection *trackcuts, double minJetPt, bool isMC, bool usePatches, bool isSwapEta){
  EMCalTriggerPtAnalysis::AliEMCalTriggerRecJetAnalysisComponent *jetana = new EMCalTriggerPtAnalysis::AliEMCalTriggerRecJetAnalysisComponent(Form("RecJetAna%f", minJetPt));
  jetana->SetMinimumJetPt(minJetPt);
  jetana->SetUsePatches(usePatches);
  jetana->SetSingleTrackCuts(trackcuts);
  //jetana->SetComponentDebugLevel(2);
  group->AddAnalysisComponent(jetana);
}

void CreateJetPtBinning(EMCalTriggerPtAnalysis::AliAnalysisTaskPtEMCalTriggerV1 *task){
  // Linear binning in steps of 10 GeV/c up to 200 GeV/c
  TArrayD binlimits(21);
  for(int i = 0; i < 21; i++) binlimits[i] = 10.*i;
  task->SetBinning("jetpt", binlimits);
}

EMCalTriggerPtAnalysis::AliEMCalPtTaskVTrackSelection *CreateDefaultTrackCuts(bool isAOD){
  EMCalTriggerPtAnalysis::AliEMCalPtTaskVTrackSelection * trackSelection(NULL);
  if(isAOD){
	  EMCalTriggerPtAnalysis::AliEMCalPtTaskTrackSelectionAOD *aodsel = new EMCalTriggerPtAnalysis::AliEMCalPtTaskTrackSelectionAOD();
	  aodsel->AddFilterBit(AliAODTrack::kTrkGlobal);
	  trackSelection = aodsel;
  } else {
	  AliESDtrackCuts *standardTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
	  standardTrackCuts->SetName("Standard Track cuts");
	  standardTrackCuts->SetMinNCrossedRowsTPC(120);
	  standardTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
	  trackSelection = new EMCalTriggerPtAnalysis::AliEMCalPtTaskTrackSelectionESD(standardTrackCuts);
  }
  return trackSelection;
}

EMCalTriggerPtAnalysis::AliEMCalPtTaskVTrackSelection *CreateHybridTrackCuts(bool isAOD){
  EMCalTriggerPtAnalysis::AliEMCalPtTaskVTrackSelection * trackSelection(NULL);
  if(isAOD){
	  // Purely use filter bits
	  EMCalTriggerPtAnalysis::AliEMCalPtTaskTrackSelectionAOD *aodsel = new EMCalTriggerPtAnalysis::AliEMCalPtTaskTrackSelectionAOD(NULL);
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
	  trackSelection = new EMCalTriggerPtAnalysis::AliEMCalPtTaskTrackSelectionESD(hybridTrackCuts);
  }
  return trackSelection;
}

EMCalTriggerPtAnalysis::AliEMCalPtTaskVTrackSelection *TrackCutsFactory(const char* trackCutsName, bool isAOD) {
  if(!strcmp(trackCutsName, "standard")) return CreateDefaultTrackCuts(isAOD);
  else if(!strcmp(trackCutsName, "hybrid")) return CreateHybridTrackCuts(isAOD);
  return NULL;
}
