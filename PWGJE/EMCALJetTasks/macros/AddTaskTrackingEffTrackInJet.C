#if ! defined __CINT__ || defined __MAKECINT__
#include <TList.h>
#include <TString.h>

#include "AliAnalysisTaskPtEfficiencyJets.h"
#include "AliAnalysisManager.h"
#include "AliEmcalTrackSelectionESD.h"
#include "AliESDtrackCuts.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#endif

PWGJE::EMCALJetTasks::AliAnalysisTaskPtEfficiencyJets *AddTaskTrackingEffTrackInJet(
    const char * mcjetcontainer  = "",
    double jetradius
    )
{
  PWGJE::EMCALJetTasks::AliAnalysisTaskPtEfficiencyJets *analysis =
      new PWGJE::EMCALJetTasks::AliAnalysisTaskPtEfficiencyJets("HighPtTrackingEfficiencyTask");

  analysis->SetIsPythia(true);

  // Connect MC jet container
  const char *jetcontainername = "HighPtEfficiencyTaskMCJetContainer";
  AliParticleContainer *mcpartcont = analysis->AddParticleContainer("MCParticlesSelected");
  AliJetContainer *jetcontainerMC = analysis->AddJetContainer(mcjetcontainer, "TPC", jetradius);
  analysis->SetMCJetContainer("jetcontainername");
  jetcontainerMC->SetName("jetcontainername");
  jetcontainerMC->ConnectParticleContainer(mcpartcont);
  jetcontainerMC->SetJetPtCut(20.);
  jetcontainerMC->SetJetEtaLimits(-1.5, 1.5);       // For this purpose the jet doesn't need to be fully in the acceptance

  // Create standard track cuts
  AliESDtrackCuts *standardTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
  standardTrackCuts->SetName("Standard Track cuts");
  standardTrackCuts->SetMinNCrossedRowsTPC(120);
  standardTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  analysis->SetTrackCuts(new AliEmcalTrackSelectionESD(standardTrackCuts));

  // connect containers
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString containerName = mgr->GetCommonFileName();
  containerName += ":HighPtTrackingEfficiency";
  printf("container name: %s\n", containerName.Data());
  mgr->ConnectInput(analysis, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(analysis, 1, mgr->CreateContainer("results", TList::Class(),    AliAnalysisManager::kOutputContainer, containerName.Data()));


  return analysis;
}
