#include <iostream>
#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGeoGlobalMagField.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliAODEvent.h"
#include "AliHeader.h"
#include "AliMCEvent.h"
#include "AliGenEventHeader.h"
#include "AliESDtrackCuts.h"
#include "AlidNdPtTools.h"
#include "AliAnalysisTaskMKBase.h"
#include "AliAnalysisTaskCutStudies.h"

class AliAnalysisTaskCutStudies;

using namespace std;

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
  fOutputList->Add(myHist.GenerateHist("myHist"));
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
 * Analyze the event.
 */
//****************************************************************************************
void AliAnalysisTaskCutStudies::AnaEvent()
{
   LoopOverAllTracks();
   if (fIsMC) LoopOverAllParticles();


}

//****************************************************************************************
/**
 * Analyze the track..
 */
//****************************************************************************************
void AliAnalysisTaskCutStudies::AnaTrack(Int_t flag)
{
    if (!fAcceptTrackM) return;

  myHist.Fill({1.,2.});
}

//****************************************************************************************
/**
 * Analyze the MC track.
 */
//****************************************************************************************
void AliAnalysisTaskCutStudies::AnaTrackMC(Int_t flag)
{
    if (!fAcceptTrackM) return;

}

//****************************************************************************************
/**
 * Analyze the MC particle.
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
        ::Error("AddTaskCutStudies", "This task requires an input event handler");
        return nullptr;
    }

    AliAnalysisTaskCutStudies *task = new AliAnalysisTaskCutStudies(name);
    if (!task) { return nullptr; }
    
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
//     task->SetESDtrackCuts(0,AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));

    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root"));
    
  return task;
}
