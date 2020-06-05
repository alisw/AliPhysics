#include <iostream>
#include "AlidNdPtTools.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskCutStudies.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskCutStudies)
/// \endcond

//****************************************************************************************
/**
 * Default constructor.
 */
//****************************************************************************************
AliAnalysisTaskCutStudies::AliAnalysisTaskCutStudies()
: AliAnalysisTaskMKBase(), fHist_x{}, fHist_y{}, fHist_z{}, fHist_alpha{}, fHist_signed1Pt{}, fHist_snp{}, fHist_tgl{},
fHist_dcaxy{}, fHist_dcaz{}, fHist_flag{}, fHist_pt{}, fHist_eta{}, fHist_phi{}, fHist_itsFoundClusters{},
fHist_itsChi2PerCluster{}, fHist_itsHits{}, fHist_tpcFindableClusters{}, fHist_tpcFoundClusters{},
fHist_tpcSharedClusters{}, fHist_tpcFractionSharedClusters{}, fHist_tpcCrossedRows{},
fHist_tpcCrossedRowsOverFindableClusters{}, fHist_tpcChi2PerCluster{}
{
}



//****************************************************************************************
/**
 * Named constructor.
 */
//****************************************************************************************
AliAnalysisTaskCutStudies::AliAnalysisTaskCutStudies(const char* name)
: AliAnalysisTaskMKBase(name), fHist_x{}, fHist_y{}, fHist_z{}, fHist_alpha{}, fHist_signed1Pt{}, fHist_snp{}, fHist_tgl{},
fHist_dcaxy{}, fHist_dcaz{}, fHist_flag{}, fHist_pt{}, fHist_eta{}, fHist_phi{}, fHist_itsFoundClusters{},
fHist_itsChi2PerCluster{}, fHist_itsHits{}, fHist_tpcFindableClusters{}, fHist_tpcFoundClusters{},
fHist_tpcSharedClusters{}, fHist_tpcFractionSharedClusters{}, fHist_tpcCrossedRows{},
fHist_tpcCrossedRowsOverFindableClusters{}, fHist_tpcChi2PerCluster{}
{
}

//****************************************************************************************
/**
 * Destructor.
 */
//****************************************************************************************
AliAnalysisTaskCutStudies::~AliAnalysisTaskCutStudies()
{
}

//****************************************************************************************
/**
 * Add output to this task.
 */
//****************************************************************************************
void AliAnalysisTaskCutStudies::AddOutput()
{
  //fHist_chi2PerClusterTPC.AddAxis("cutSetting", "[noCuts, defaultCut]", 2, -0.5, 1.5);
  fHist_tpcChi2PerCluster.AddAxis("chi2PerClusterTPC", "chi2/cls TPC", 500, 0., 50.);
  fOutputList->Add(fHist_tpcChi2PerCluster.GenerateHist("tpc-chi2PerCluster"));
  
  fHist_tpcFoundClusters.AddAxis("foundClustersTPC", "# found clusters TPC",  165, -0.5, 164.5);
  fOutputList->Add(fHist_tpcFoundClusters.GenerateHist("tpc-foundClusters"));

  fHist_tpcFindableClusters.AddAxis("findableClustersTPC", "# findable clusters TPC",  165, -0.5, 164.5);
  fOutputList->Add(fHist_tpcFindableClusters.GenerateHist("tpc-findableClusters"));

  fHist_tpcCrossedRows.AddAxis("crossedRowsTPC", "# crossed rows TPC",  165, -0.5, 164.5);
  fOutputList->Add(fHist_tpcCrossedRows.GenerateHist("tpc-crossedRows"));
  
  fHist_tpcCrossedRowsOverFindableClusters.AddAxis("crossedRowsOverFindableTPC", "# crossed rows per findable TPC",  165, -0.5, 164.5);
  fOutputList->Add(fHist_tpcCrossedRowsOverFindableClusters.GenerateHist("tpc-crossedRows"));
  
  
}

//****************************************************************************************
/**
 * Event selection.
 */
//****************************************************************************************
Bool_t AliAnalysisTaskCutStudies::IsEventSelected()
{
  return fIsAcceptedAliEventCuts;
}

//****************************************************************************************
/**
 * Analyse the event.
 */
//****************************************************************************************
void AliAnalysisTaskCutStudies::AnaEvent()
{
   LoopOverAllTracks();
   if (fIsMC) LoopOverAllParticles();
}

//****************************************************************************************
/**
 * Analyse the track.
 */
//****************************************************************************************
void AliAnalysisTaskCutStudies::AnaTrack(Int_t flag)
{
  
  fHist_tpcChi2PerCluster.Fill(10.);
  fHist_tpcFoundClusters.Fill(10.);
  fHist_tpcFindableClusters.Fill(10.);
  fHist_tpcCrossedRows.Fill(10.);
  fHist_tpcCrossedRowsOverFindableClusters.Fill(10.);


  //if (!fAcceptTrackM) return;
  
  /*
  for(int i=0; i<8;++i)
  {
    if(fAcceptTrack[i])
    {
      fHist_chi2PerClusterTPC.Fill(i, 10.);

    }
  }
   */
}

//****************************************************************************************
/**
 * Analyse the MC track.
 */
//****************************************************************************************
void AliAnalysisTaskCutStudies::AnaTrackMC(Int_t flag)
{
    if (!fAcceptTrackM) return;

}

//****************************************************************************************
/**
 * Analyse the MC particle.
 */
//****************************************************************************************
void AliAnalysisTaskCutStudies::AnaParticleMC(Int_t flag)
{            
    if (!fMCisPrim) return;    
    if (!fMCIsCharged) return;    
    if (TMath::Abs(fMCEta) > 0.8) return;    
  
}

//****************************************************************************************
/**
 * Add task of this kind to a train.
 */
//****************************************************************************************
AliAnalysisTaskCutStudies* AliAnalysisTaskCutStudies::AddTaskCutStudies(const char* name)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
      ::Error("AddTaskCutStudies", "No analysis manager to connect to.");
      return nullptr;
  }

  if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskCutStudies", "This task requires an input event handler.");
      return nullptr;
  }

  AliAnalysisTaskCutStudies *task = new AliAnalysisTaskCutStudies(name);
  if (!task) { return nullptr; }

  task->SelectCollisionCandidates(AliVEvent::kAnyINT);
  task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root"));

  return task;
}
