#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskMKBase.h"
#include "AliAnalysisTaskSpectraEtaPhi.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliGenEventHeader.h"
#include "AliHeader.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AlidNdPtTools.h"
#include "TChain.h"
#include "TGeoGlobalMagField.h"
#include "TH1F.h"
#include "TList.h"
#include "TRandom.h"
#include "TRandom3.h"
#include <iostream>

namespace {
using namespace AnalysisHelpers;
using namespace std;
} // namespace

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSpectraEtaPhi)
    /// \endcond
    //_____________________________________________________________________________

    AliAnalysisTaskSpectraEtaPhi::AliAnalysisTaskSpectraEtaPhi()
    : AliAnalysisTaskMKBase(), fHistEffContNCluster{}, fHistEffContZ{},
      fHistEffContEta{}, fHistEffContPhi{}, fHistTrackNCluster{}, fHistTrackZ{},
      fHistTrackEta{}, fHistTrackPhi{}, fHistEvent{} {
    // default contructor
}

//_____________________________________________________________________________

AliAnalysisTaskSpectraEtaPhi::AliAnalysisTaskSpectraEtaPhi(const char* name)
    : AliAnalysisTaskMKBase(name), fHistEffContNCluster{}, fHistEffContZ{},
      fHistEffContEta{}, fHistEffContPhi{}, fHistTrackNCluster{}, fHistTrackZ{},
      fHistTrackEta{}, fHistTrackPhi{}, fHistEvent{} {
    // constructor
}

//_____________________________________________________________________________

AliAnalysisTaskSpectraEtaPhi::~AliAnalysisTaskSpectraEtaPhi() {
    // destructor
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraEtaPhi::AddOutput() {
    // general binnings
    std::vector<double> centBins = {0., 10, 20., 30., 40, 50, 60, 70, 80, 90};
    std::vector<double> ptBins = {0.0, 0.1, 0.2, 0.3, 0.4,  0.5, 0.6, 0.7,
                                  0.8, 0.9, 1.0, 1.1, 1.2,  1.3, 1.4, 1.5,
                                  2.0, 3.5, 5.0, 7.5, 10.0, 20.0};
    const Int_t nbins = 54;
    std::vector<double> multBins;
    multBins.reserve(nbins + 1);
    multBins.push_back(-0.5);
    {
        int i = 0;
        for (; i <= 10; i++) {
            multBins.push_back(multBins.back() + 1);
        }
        for (; i <= 10 + 9; i++) {
            multBins.push_back(multBins.back() + 10);
        }
        for (; i <= 10 + 9 + 9; i++) {
            multBins.push_back(multBins.back() + 100);
        }
        for (; i <= 10 + 9 + 9 + 25; i++) {
            multBins.push_back(multBins.back() + 200);
        }
    }

    const int nCuts = 8;
    Axis cutAxis = {"cut", "cut setting", {-1.5, nCuts - 0.5}, nCuts + 1};
    Axis centAxis = {"cent", "centrality", centBins};
    Axis ptAxis = {"pt", "#it{p}_{T} (GeV/c)", ptBins};
    Axis etaAxis = {"eta", "#eta", {-0.8, 0.8}, 8};
    Axis phiAxis = {
        "phi", "#phi", {0., 2. * M_PI}, 36}; // 36 to see tpc sectors
    Axis mcInfoAxis = {"mcInfo", "mcInfo", {-0.3, 3.5}, 4};
    Axis NClusterAxis = {"NClusterPID", "N_{Cluster}^{PID}", {-0.5, 160.5}, 20};
    Axis ZInnerAxis = {"ZInnerParam", "Z_{Inner}", {-30.5, 29.5}, 60};
    Axis ChargeQAxis = {"MCQ", "Q", {-1.5, 1.5}, 3};
    Axis multAxis = {"mult", "#it{N}_{ch}", multBins};
    Axis zVrtAxis = {"zV", "vertex_{Z}", {-20, 20}, 8};

    double requiredMemory = 0;
    // MC
    fHistEffContNCluster.AddAxis(centAxis);
    fHistEffContNCluster.AddAxis(ptAxis);
    fHistEffContNCluster.AddAxis(NClusterAxis);
    fHistEffContNCluster.AddAxis(cutAxis);
    fHistEffContNCluster.AddAxis(mcInfoAxis);
    fOutputList->Add(fHistEffContNCluster.GenerateHist("fHistEffContNCluster"));
    requiredMemory += fHistEffContNCluster.GetSize();

    fHistEffContZ.AddAxis(centAxis);
    fHistEffContZ.AddAxis(ptAxis);
    fHistEffContZ.AddAxis(ZInnerAxis);
    fHistEffContZ.AddAxis(cutAxis);
    fHistEffContZ.AddAxis(mcInfoAxis);
    fOutputList->Add(fHistEffContZ.GenerateHist("fHistEffContZ"));
    requiredMemory += fHistEffContZ.GetSize();

    fHistEffContEta.AddAxis(centAxis);
    fHistEffContEta.AddAxis(ptAxis);
    fHistEffContEta.AddAxis(etaAxis);
    fHistEffContEta.AddAxis(cutAxis);
    fHistEffContEta.AddAxis(mcInfoAxis);
    fHistEffContEta.AddAxis(ChargeQAxis);
    fOutputList->Add(fHistEffContEta.GenerateHist("fHistEffContEta"));
    requiredMemory += fHistEffContEta.GetSize();

    fHistEffContPhi.AddAxis(centAxis);
    fHistEffContPhi.AddAxis(ptAxis);
    fHistEffContPhi.AddAxis(phiAxis);
    fHistEffContPhi.AddAxis(cutAxis);
    fHistEffContPhi.AddAxis(mcInfoAxis);
    fOutputList->Add(fHistEffContPhi.GenerateHist("fHistEffContPhi"));
    requiredMemory += fHistEffContPhi.GetSize();

    // data
    fHistTrackNCluster.AddAxis(centAxis);
    fHistTrackNCluster.AddAxis(ptAxis);
    fHistTrackNCluster.AddAxis(NClusterAxis);
    fHistTrackNCluster.AddAxis(cutAxis);
    fOutputList->Add(fHistTrackNCluster.GenerateHist("fHistTrackNCluster"));
    requiredMemory += fHistTrackNCluster.GetSize();

    fHistTrackZ.AddAxis(centAxis);
    fHistTrackZ.AddAxis(ptAxis);
    fHistTrackZ.AddAxis(ZInnerAxis);
    fHistTrackZ.AddAxis(cutAxis);
    fOutputList->Add(fHistTrackZ.GenerateHist("fHistTrackZ"));
    requiredMemory += fHistTrackZ.GetSize();

    fHistTrackEta.AddAxis(centAxis);
    fHistTrackEta.AddAxis(ptAxis);
    fHistTrackEta.AddAxis(etaAxis);
    fHistTrackEta.AddAxis(cutAxis);
    fHistTrackEta.AddAxis(ChargeQAxis);
    fOutputList->Add(fHistTrackEta.GenerateHist("fHistTrackEta"));
    requiredMemory += fHistTrackEta.GetSize();

    fHistTrackPhi.AddAxis(centAxis);
    fHistTrackPhi.AddAxis(ptAxis);
    fHistTrackPhi.AddAxis(phiAxis);
    fHistTrackPhi.AddAxis(cutAxis);
    fOutputList->Add(fHistTrackPhi.GenerateHist("fHistTrackPhi"));
    requiredMemory += fHistTrackPhi.GetSize();

    // Event hists
    fHistEvent.AddAxis(centAxis);
    fHistEvent.AddAxis(multAxis);
    fHistEvent.AddAxis(zVrtAxis);
    fOutputList->Add(fHistEvent.GenerateHist("fHistEvent"));
    requiredMemory += fHistEvent.GetSize();

    AliError(Form("Estimated memory usage of histograms: %.0f Bytes (%f MiB)",
                  requiredMemory, requiredMemory / 1048576));
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskSpectraEtaPhi::IsEventSelected() {
    return fIsAcceptedAliEventCuts;
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraEtaPhi::AnaEvent() {

    LoopOverAllTracks();
    if (fIsMC)
        LoopOverAllParticles();

    fHistEvent.Fill(fMultPercentileV0M, fNTracksAcc, fZv);
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraEtaPhi::AnaTrack(Int_t flag) {
    //    if (!fAcceptTrackM) return;
    for (int i = 0; i < 8; ++i) {
        if (fAcceptTrack[i]) {
            fHistTrackNCluster.Fill(fMultPercentileV0M, fPt, fTPCSignalN, i);
            fHistTrackZ.Fill(fMultPercentileV0M, fPt, fZInner, i);
            fHistTrackEta.Fill(fMultPercentileV0M, fPt, fEta, i, fChargeSign);
            fHistTrackPhi.Fill(fMultPercentileV0M, fPt, fPhi, i);
        }
    }
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraEtaPhi::AnaTrackMC(Int_t flag) {
    //    if (!fAcceptTrackM) return;

    if (fMCParticleType == AlidNdPtTools::kOther) {
        Log("RecTrack.PDG.", fMCPDGCode);
    }
    if (TMath::Abs(fMCQ > 1)) {
        Log("RecTrack.Q>1.PDG.", fMCPDGCode);
    }

    for (int i = 0; i < 8; ++i) {
        if (fAcceptTrack[i]) {
            fHistEffContNCluster.Fill(fMultPercentileV0M, fPt, fTPCSignalN, i,
                                      fMCProdcutionType);
            fHistEffContZ.Fill(fMultPercentileV0M, fPt, fZInner, i,
                               fMCProdcutionType);
            fHistEffContEta.Fill(fMultPercentileV0M, fPt, fEta, i,
                                 fMCProdcutionType, fChargeSign);
            fHistEffContPhi.Fill(fMultPercentileV0M, fPt, fPhi, i,
                                 fMCProdcutionType);
        }
    }
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraEtaPhi::AnaParticleMC(Int_t flag) {
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

    fHistEffContNCluster.Fill(fMultPercentileV0M, fMCPt, fTPCSignalN, -1., 3.);
    fHistEffContZ.Fill(fMultPercentileV0M, fMCPt, fZInner, -1., 3.);
    fHistEffContEta.Fill(fMultPercentileV0M, fMCPt, fMCEta, -1., 3., fChargeSign);
    fHistEffContPhi.Fill(fMultPercentileV0M, fMCPt, fMCPhi, -1., 3.);
}

//_____________________________________________________________________________

AliAnalysisTaskSpectraEtaPhi*
AliAnalysisTaskSpectraEtaPhi::AddTaskSpectra(const char* name,
                                             const char* outfile) {
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskSpectraEtaPhi", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the
    // analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskSpectraEtaPhi",
                "This task requires an input event handler");
        return NULL;
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
    AliAnalysisTaskSpectraEtaPhi* task = new AliAnalysisTaskSpectraEtaPhi(name);
    if (!task) {
        return 0;
    }

    // configure the task
    //===========================================================================
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
    // default, nogeo, nogeonogold, geoCut116, geoCut117, geoCut118 geoCut119,
    // tpconlyminimal
    task->SetESDtrackCuts(0, AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
    task->SetESDtrackCuts(
        1, AlidNdPtTools::CreateESDtrackCuts("tpcitsnogeoEta08"));
    task->SetESDtrackCuts(
        2, AlidNdPtTools::CreateESDtrackCuts("tpcitsnogeonogoldEta08"));
    task->SetESDtrackCuts(
        3, AlidNdPtTools::CreateESDtrackCuts("defaultEta08", 116));
    task->SetESDtrackCuts(
        4, AlidNdPtTools::CreateESDtrackCuts("defaultEta08", 117));
    task->SetESDtrackCuts(
        5, AlidNdPtTools::CreateESDtrackCuts("defaultEta08", 118));
    task->SetESDtrackCuts(
        6, AlidNdPtTools::CreateESDtrackCuts("defaultEta08", 119));
    task->SetESDtrackCuts(
        7, AlidNdPtTools::CreateESDtrackCuts("tpconlyminimalEta08"));
    task->SetNeedEventMult(kTRUE);
    task->SetNeedTrackIP(kTRUE);

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
