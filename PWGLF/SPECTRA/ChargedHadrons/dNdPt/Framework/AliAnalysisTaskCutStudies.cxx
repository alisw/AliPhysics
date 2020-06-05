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
: AliAnalysisTaskMKBase(), myHist{}
{
}

//****************************************************************************************
/**
 * Named constructor.
 */
//****************************************************************************************
AliAnalysisTaskCutStudies::AliAnalysisTaskCutStudies(const char* name)
: AliAnalysisTaskMKBase(name), myHist{}
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
  myHist.AddAxis("axis1", "title1", 100, 0, 100);
  myHist.AddAxis("axis2", "title2", 100, 0, 100);
  fOutputList->Add(myHist.GenerateHist("myHistABC"));
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
  if (!fAcceptTrackM) return;
  myHist.Fill(10., 20.);
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
