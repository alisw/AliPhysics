#if !defined (__CINT__) || defined (__MAKECINT__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskPtEMCalTriggerV1.h"
#include "AliESDtrackCuts.h"
#include "AliJetContainer.h"
#include <TList.h>
#include <TString.h>
#include <cstring>
#endif

void AddClusterComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, bool isMC);
void AddTrackComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, AliESDtrackCuts *trackcuts, bool isMC, bool isSwapEta);
void AddMCParticleComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group);
void AddEventCounterComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, bool isMC);
void AddMCJetComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, double minJetPt);
void AddRecJetComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, AliESDtrackCuts *trackcuts, double minJetPt, bool isMC, bool isSwapEta);
AliESDtrackCuts *CreateDefaultTrackCuts();
AliESDtrackCuts *CreateHybridTrackCuts();

AliAnalysisTask* AddTaskPtEMCalTriggerV1(
    bool isMC,
    bool usePythiaHard,
    const char *period ="LHC13d",
    const char *ntrackContainer = "",
    const char *nclusterContainer = "",
    const char *njetcontainerData = "",
    const char *njetcontainerMC = "",
    double jetradius = 0.5
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

  bool isSwapEta = TString(period).CompareTo("LHC13f") ? kFALSE : kTRUE;
  EMCalTriggerPtAnalysis::AliAnalysisTaskPtEMCalTriggerV1 *pttriggertask = new EMCalTriggerPtAnalysis::AliAnalysisTaskPtEMCalTriggerV1("ptemcaltriggertask");
  //pttriggertask->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kEMC7);                          // Select both INT7 or EMC7 triggered events
  pttriggertask->SelectCollisionCandidates(AliVEvent::kAny);
  if(isMC) pttriggertask->SetSwapThresholds();

  mgr->AddTask(pttriggertask);
  if(usePythiaHard){
    pttriggertask->SetIsPythia(kTRUE);
  }

  // Add components
  EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *noselect = new EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup("noselect");
  noselect->AddAnalysisComponent(new EMCalTriggerPtAnalysis::AliEMCalTriggerPatchAnalysisComponent("patchanalysis"));

  double jetpt[4] = {40., 60., 80., 100.};
  EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *defaultselect = new EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup("defaultselect");
  defaultselect->SetEventSelection(new EMCalTriggerPtAnalysis::AliEMCalTriggerEventSelection());
  EMCalTriggerPtAnalysis::AliEMCalTriggerKineCuts *kineCuts = new EMCalTriggerPtAnalysis::AliEMCalTriggerKineCuts();
  kineCuts->SetPtRange(2., 100.);
  defaultselect->SetKineCuts(kineCuts);
  AddEventCounterComponent(defaultselect, isMC);
  if(isMC){
    AddMCParticleComponent(defaultselect);
    for(int ijpt = 0; ijpt < 4; ijpt++)
      AddMCJetComponent(defaultselect, jetpt[ijpt]);
  }
  AddClusterComponent(defaultselect, isMC);
  AddTrackComponent(defaultselect, CreateDefaultTrackCuts(), isMC, isSwapEta);
  for(int ijpt = 0; ijpt < 4; ijpt++)
    AddRecJetComponent(defaultselect, CreateDefaultTrackCuts(), jetpt[ijpt], isMC, isSwapEta);

  pttriggertask->AddAnalysisGroup(noselect);
  pttriggertask->AddAnalysisGroup(defaultselect);

  // Add containers
  Bool_t isAOD = mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class();
  AliParticleContainer *trackContainer = pttriggertask->AddParticleContainer(ntrackContainer);
  //trackContainer->SetClassName("AliVTrack");
  AliClusterContainer *clusterContainer = pttriggertask->AddClusterContainer(nclusterContainer);
  AliParticleContainer *mcpartcont = isMC ? pttriggertask->AddParticleContainer("MCParticlesSelected") : NULL;


  // Handle Jet Containers
  if(strlen(njetcontainerData)){
    AliJetContainer *jetcontainerData = pttriggertask->AddJetContainer(njetcontainerData, "TPC", jetradius);
    pttriggertask->SetDataJetContainerName("PtTriggerTaskJetsData");
    jetcontainerData->ConnectParticleContainer(trackContainer);
    jetcontainerData->SetName("PtTriggerTaskJetsData");
    jetcontainerData->SetJetPtCut(20.);
    jetcontainerData->SetJetEtaPhiEMCAL();
    printf("jet container added for Data\n");
  }
  if(isMC && strlen(njetcontainerMC)){
    AliJetContainer *jetcontainerMC = pttriggertask->AddJetContainer(njetcontainerMC, "TPC", jetradius);
    pttriggertask->SetMCJetContainerName("PtTriggerTaskJetsMC");
    jetcontainerMC->ConnectParticleContainer(mcpartcont);
    jetcontainerMC->SetName("PtTriggerTaskJetsMC");
    jetcontainerMC->SetJetPtCut(20.);
    jetcontainerMC->SetJetEtaPhiEMCAL();
    printf("Jet container added for MC");
  }

  TString containerName = mgr->GetCommonFileName();
  containerName += ":PtEMCalTriggerTask";
  printf("container name: %s\n", containerName.Data());

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("results", TList::Class(),    AliAnalysisManager::kOutputContainer, containerName.Data());

  //Connect input/output
  mgr->ConnectInput(pttriggertask, 0, cinput);
  mgr->ConnectOutput(pttriggertask, 1, coutput);

  return pttriggertask;
}

void AddClusterComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, bool isMC){
  EMCalTriggerPtAnalysis::AliEMCalTriggerClusterAnalysisComponent *clusteranalysis = new EMCalTriggerPtAnalysis::AliEMCalTriggerClusterAnalysisComponent("clusterAnalysis");
  clusteranalysis->SetEnergyRange(2., 100.);
  if(isMC) clusteranalysis->SetUsePatches();
  group->AddAnalysisComponent(clusteranalysis);
}

void AddTrackComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, AliESDtrackCuts * trackcuts, bool isMC, bool isSwapEta){
  EMCalTriggerPtAnalysis::AliEMCalTriggerRecTrackAnalysisComponent *trackanalysis = new EMCalTriggerPtAnalysis::AliEMCalTriggerRecTrackAnalysisComponent("trackAnalysisStandard");
  group->AddAnalysisComponent(trackanalysis);
  // Create charged hadrons pPb standard track cuts
  trackanalysis->SetTrackSelection(new EMCalTriggerPtAnalysis::AliEMCalPtTaskTrackSelectionESD(trackcuts));

  if(isMC){
    trackanalysis->SetRequestMCtrueTracks();
    trackanalysis->SetUsePatches();
  }
  if(isSwapEta) trackanalysis->SetSwapEta();
}

void AddEventCounterComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, Bool_t isMC){
  EMCalTriggerPtAnalysis::AliEMCalTriggerEventCounterAnalysisComponent * evcount = new EMCalTriggerPtAnalysis::AliEMCalTriggerEventCounterAnalysisComponent("eventCounter");
  evcount->SetUsePatches(isMC);
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

void AddRecJetComponent(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup *group, AliESDtrackCuts *trackcuts, double minJetPt, bool isMC, bool isSwapEta){
  EMCalTriggerPtAnalysis::AliEMCalTriggerRecJetAnalysisComponent *jetana = new EMCalTriggerPtAnalysis::AliEMCalTriggerRecJetAnalysisComponent(Form("RecJetAna%f", minJetPt));
  jetana->SetMinimumJetPt(minJetPt);
  jetana->SetUsePatches();
  jetana->SetSingleTrackCuts(new EMCalTriggerPtAnalysis::AliEMCalPtTaskTrackSelectionESD(trackcuts));
  group->AddAnalysisComponent(jetana);
}

AliESDtrackCuts *CreateDefaultTrackCuts(){
  AliESDtrackCuts *standardTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
  standardTrackCuts->SetName("Standard Track cuts");
  standardTrackCuts->SetMinNCrossedRowsTPC(120);
  standardTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  return standardTrackCuts;
}

AliESDtrackCuts *CreateHybridTrackCuts(){
  AliESDtrackCuts* hybridTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
  hybridTrackCuts->SetName("Global Hybrid tracks, loose DCA");
  hybridTrackCuts->SetMaxDCAToVertexXY(2.4);
  hybridTrackCuts->SetMaxDCAToVertexZ(3.2);
  hybridTrackCuts->SetDCAToVertex2D(kTRUE);
  hybridTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
  hybridTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
  return hybridTrackCuts;
}
