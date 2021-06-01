#if !defined (__CINT__) || defined (__MAKECINT__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskPtEMCalTrigger.h"
#include "AliESDtrackCuts.h"
#include "AliJetContainer.h"
#include <TList.h>
#include <TString.h>
#include <cstring>
#endif

AliAnalysisTask* AddTaskPtEMCalTrigger(
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

  PWGJE::EMCALJetTasks::AliAnalysisTaskPtEMCalTrigger *pttriggertask = new PWGJE::EMCALJetTasks::AliAnalysisTaskPtEMCalTrigger("ptemcaltriggertask");
  //pttriggertask->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kEMC7);                          // Select both INT7 or EMC7 triggered events
  pttriggertask->SelectCollisionCandidates(AliVEvent::kAny);
  if(!TString(period).CompareTo("LHC13f")) pttriggertask->SetSwapEta();
  mgr->AddTask(pttriggertask);
  pttriggertask->SetPtRange(2., 100.);
  pttriggertask->SetClusterEnergyRange(2.,100.);
  if(usePythiaHard){
    pttriggertask->SetIsPythia(kTRUE);
  }

  // Add containers
  Bool_t isAOD = mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class();
  AliParticleContainer *trackContainer = pttriggertask->AddParticleContainer(ntrackContainer);
  //trackContainer->SetClassName("AliVTrack");
  AliClusterContainer *clusterContainer = pttriggertask->AddClusterContainer(nclusterContainer);
  AliParticleContainer *mcpartcont = isMC ? pttriggertask->AddParticleContainer("MCParticlesSelected") : NULL;


  // Create charged hadrons pPb standard track cuts
  AliESDtrackCuts *standardTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
  standardTrackCuts->SetName("Standard Track cuts");
  standardTrackCuts->SetMinNCrossedRowsTPC(120);
  standardTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  pttriggertask->AddESDTrackCuts(standardTrackCuts);

  // Create hybrid track cuts as used in the jet analysis
  AliESDtrackCuts* hybridTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
  hybridTrackCuts->SetName("Global Hybrid tracks, loose DCA");
  hybridTrackCuts->SetMaxDCAToVertexXY(2.4);
  hybridTrackCuts->SetMaxDCAToVertexZ(3.2);
  hybridTrackCuts->SetDCAToVertex2D(kTRUE);
  hybridTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
  hybridTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
  pttriggertask->AddESDTrackCuts(hybridTrackCuts);

  // Handle Jet Containers
  if(strlen(njetcontainerData)){
    AliJetContainer *jetcontainerData = pttriggertask->AddJetContainer(njetcontainerData, "TPC", jetradius);
    pttriggertask->AddJetContainerName("PtTriggerTaskJetsData", false);
    jetcontainerData->ConnectParticleContainer(trackContainer);
    jetcontainerData->SetName("PtTriggerTaskJetsData");
    jetcontainerData->SetJetPtCut(20.);
  }
  if(isMC && strlen(njetcontainerMC)){
    AliJetContainer *jetcontainerMC = pttriggertask->AddJetContainer(njetcontainerMC, "TPC", jetradius);
    pttriggertask->AddJetContainerName("PtTriggerTaskJetsMC", true);
    jetcontainerMC->ConnectParticleContainer(mcpartcont);
    jetcontainerMC->SetName("PtTriggerTaskJetsMC");
    jetcontainerMC->SetJetPtCut(20.);
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
