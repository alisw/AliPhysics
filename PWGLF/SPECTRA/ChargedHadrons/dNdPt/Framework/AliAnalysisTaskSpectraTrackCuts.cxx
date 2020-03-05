#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskMKBase.h"
#include "AliAnalysisTaskSpectraTrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliGenEventHeader.h"
#include "AliHeader.h"
#include "AliMCEvent.h"
#include "AlidNdPtTools.h"
#include "TChain.h"
#include "TGeoGlobalMagField.h"
#include "TH1F.h"
#include "TList.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TFile.h"
#include <iostream>

class AliAnalysisTaskSpectraTrackCuts;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSpectraTrackCuts)
    /// \endcond
    //_____________________________________________________________________________

    AliAnalysisTaskSpectraTrackCuts::AliAnalysisTaskSpectraTrackCuts()
    : AliAnalysisTaskMKBase() {
    // default contructor
}

//_____________________________________________________________________________

AliAnalysisTaskSpectraTrackCuts::AliAnalysisTaskSpectraTrackCuts(
    const char* name)
    : AliAnalysisTaskMKBase(name) {
    // constructor
}

//_____________________________________________________________________________

AliAnalysisTaskSpectraTrackCuts::~AliAnalysisTaskSpectraTrackCuts() {
    // destructor
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraTrackCuts::Terminate(Option_t *option) {
    // Terminate
    std::cout << "Saving ESDTrackCuts settings\n";
    gFile->cd();
    std::cout << "Currently in file " << gFile->GetName() << "\n";
    if (fESDtrackCutsM) {
        gFile->mkdir("ESD_TrackCuts");
        fESDtrackCutsM->SaveHistograms("ESD_TrackCuts");
    }
    for (int i = 0; i < 10; ++i) {
        if (fESDtrackCuts[i]) {
            gFile->mkdir(Form("ESD_TrackCuts_%d", i));
            fESDtrackCuts[i]->SaveHistograms(Form("ESD_TrackCuts_%d", i));
        }
    }
    std::cout << "...successful\n";
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskSpectraTrackCuts::IsEventSelected() {
    return fIsAcceptedAliEventCuts;
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraTrackCuts::AnaEvent() {

    LoopOverAllTracks();
    if (fIsMC)
        LoopOverAllParticles();
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraTrackCuts::AnaTrack(Int_t flag) {
    if (!fAcceptTrackM)
        return;
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraTrackCuts::AnaTrackMC(Int_t flag) {
    if (!fAcceptTrackM)
        return;

    if (fMCParticleType == AlidNdPtTools::kOther) {
        Log("RecTrack.PDG.", fMCPDGCode);
    }
    if (TMath::Abs(fMCQ > 1)) {
        Log("RecTrack.Q>1.PDG.", fMCPDGCode);
    }
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraTrackCuts::AnaParticleMC(Int_t flag) {
    if (!fMCisPrim)
        return;
    if (!fMCIsCharged)
        return;
    if (TMath::Abs(fMCEta) > 0.8)
        return;

    if (fMCParticleType == AlidNdPtTools::kOther) {
        Log("GenPrim.PDG.", fMCPDGCode);
    }
    if (TMath::Abs(fMCQ > 1)) {
        Log("GenPrim.Q>1.PDG.", fMCPDGCode);
    }
}

//_____________________________________________________________________________

AliAnalysisTaskSpectraTrackCuts*
AliAnalysisTaskSpectraTrackCuts::AddTaskSpectraTrackCuts(const char* name,
                                                         const char* outfile,
                                                         int _CutMode) {
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskSpectraTrackCuts",
                "No analysis manager to connect to.");
        return nullptr;
    }

    // Check the analysis type using the event handlers connected to the
    // analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskSpectraTrackCuts",
                "This task requires an input event handler");
        return nullptr;
    }

    // Setup output file
    //===========================================================================
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":";
    fileName += name; // create a subfolder in the file
    if (outfile) {    // if a finename is given, use that one
        fileName = TString(outfile);
    }

    // create the task
    //===========================================================================
    AliAnalysisTaskSpectraTrackCuts* task =
        new AliAnalysisTaskSpectraTrackCuts(name);
    if (!task) {
        return nullptr;
    }

    // configure the task
    //===========================================================================
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);
    auto const cuts =
        AlidNdPtTools::CreateESDtrackCuts("defaultEta08", _CutMode, true);
    //    cuts->SaveHistograms(fileName.Data());
    task->SetESDtrackCutsM(cuts);
    //     task->SetESDtrackCuts(0,AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));

    // attach the task to the manager and configure in and ouput
    //===========================================================================
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(
        task, 1,
        mgr->CreateContainer(name, TList::Class(),
                             AliAnalysisManager::kOutputContainer,
                             fileName.Data()));

    return task;
}
