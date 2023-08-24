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
#include "AliMCSpectraWeights.h"
#include <iostream>

using namespace std;

//_____________________________________________________________________________

AliAnalysisTaskSpectraEtaPhi::AliAnalysisTaskSpectraEtaPhi()
: AliAnalysisTaskMKBase(), fHistEffContNCluster{}, fHistEffContZ{},
  fHistEffContEta{}, fHistEffContPhi{}, fHistTrackNCluster{}, fHistTrackZ{},
  fHistTrackEta{}, fHistTrackPhi{}, fHistEvent{} {
    // default constructor
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
    std::vector<double> centBins = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
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
            multBins.push_back(multBins.back() + 1.);
        }
        for (; i <= 10 + 9; i++) {
            multBins.push_back(multBins.back() + 10.);
        }
        for (; i <= 10 + 9 + 9; i++) {
            multBins.push_back(multBins.back() + 100.);
        }
        for (; i <= 10 + 9 + 9 + 25; i++) {
            multBins.push_back(multBins.back() + 200.);
        }
    }

    const int nCuts = 23;
    Axis cutAxis = {"cut", "cut setting", {-1.5, nCuts - 0.5}, nCuts + 1};
    Axis centAxis = {"cent", "centrality", centBins};
    Axis ptAxis = {"pt", "#it{p}_{T} (GeV/c)", ptBins};
    Axis etaAxis = {"eta", "#eta", {-1, 1}, 20};
    Axis phiAxis = {
        "phi", "#phi", {0., 2. * M_PI}, 36}; // 36 to see tpc sectors
    Axis mcInfoAxis = {"mcInfo", "mcInfo", {-0.3, 3.5}, 4};
    Axis nClusterAxis = {"NClusterPID", "N_{Cluster}^{PID}", {-0.5, 160.5}, 20};
    Axis zInnerAxis = {"ZInnerParam", "Z_{Inner}", {-90.5, 90.5}, 180};
    Axis chargeQAxis = {"MCQ", "Q", {-1.5, 1.5}, 3};
    Axis multAxis = {"mult", "#it{N}_{ch}", multBins};
    Axis zVrtAxis = {"zV", "vertex_{Z}", {-20, 20}, 8};
    Axis WeightSysAxis = {"WeightSys", "WeightSys", {-1.5, 1.5}, 3};

    double requiredMemory = 0.;
    // MC
    fHistEffContNCluster.AddAxis(centAxis);
    fHistEffContNCluster.AddAxis(ptAxis);
    fHistEffContNCluster.AddAxis(nClusterAxis);
    fHistEffContNCluster.AddAxis(cutAxis);
    fHistEffContNCluster.AddAxis(mcInfoAxis);
    fOutputList->Add(fHistEffContNCluster.GenerateHist("fHistEffContNCluster"));
    requiredMemory += fHistEffContNCluster.GetSize();
    
    fHistEffContZ.AddAxis(centAxis);
    fHistEffContZ.AddAxis(ptAxis);
    fHistEffContZ.AddAxis(zInnerAxis);
    fHistEffContZ.AddAxis(cutAxis);
    fHistEffContZ.AddAxis(mcInfoAxis);
    fHistEffContZ.AddAxis(etaAxis);
    fOutputList->Add(fHistEffContZ.GenerateHist("fHistEffContZ"));
    requiredMemory += fHistEffContZ.GetSize();

    fHistEffContEta.AddAxis(multAxis);
    fHistEffContEta.AddAxis(ptAxis);
    fHistEffContEta.AddAxis(etaAxis);
    fHistEffContEta.AddAxis(cutAxis);
    fHistEffContEta.AddAxis(mcInfoAxis);
    fHistEffContEta.AddAxis(chargeQAxis);
    fHistEffContEta.AddAxis(WeightSysAxis);
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
    fHistTrackNCluster.AddAxis(nClusterAxis);
    fHistTrackNCluster.AddAxis(cutAxis);
    fOutputList->Add(fHistTrackNCluster.GenerateHist("fHistTrackNCluster"));
    requiredMemory += fHistTrackNCluster.GetSize();

    fHistTrackZ.AddAxis(centAxis);
    fHistTrackZ.AddAxis(ptAxis);
    fHistTrackZ.AddAxis(zInnerAxis);
    fHistTrackZ.AddAxis(cutAxis);
    fHistTrackZ.AddAxis(etaAxis);
    fOutputList->Add(fHistTrackZ.GenerateHist("fHistTrackZ"));
    requiredMemory += fHistTrackZ.GetSize();

    fHistTrackEta.AddAxis(centAxis);
    fHistTrackEta.AddAxis(ptAxis);
    fHistTrackEta.AddAxis(etaAxis);
    fHistTrackEta.AddAxis(cutAxis);
    fHistTrackEta.AddAxis(chargeQAxis);
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

void AliAnalysisTaskSpectraEtaPhi::AnaEventMC() {

    AliMCSpectraWeightsHandler* mcWeightsHandler = static_cast<AliMCSpectraWeightsHandler*>(fEvent->FindListObject("fMCSpectraWeights"));
    fMCSpectraWeights = (mcWeightsHandler) ? mcWeightsHandler->fMCSpectraWeight : nullptr;

    LoopOverAllParticles();
    LoopOverAllTracks();
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraEtaPhi::AnaTrack(Int_t flag) {
    //    if (!fAcceptTrackM) return;
    for (int i = 0; i < 23; ++i) {
        if (fAcceptTrack[i]) {
            fHistTrackNCluster.Fill(fMultPercentileV0M, fPt, fTPCSignalN, i);
            fHistTrackZ.Fill(fMultPercentileV0M, fPt, fZInner, i, fEta);
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
    if(fMCPileUpTrack){ // reject pileup tracks
        return;
    }

    double fMCweight = 1.0;
    double fMCweightSysUp = 1.0;
    double fMCweightSysDown = 1.0;
    if(fMCSpectraWeights && 0==fMCPrimSec && !fMCPileUpTrack && fMCParticle->Particle()){ // only for primary particles
        fMCweight = fMCSpectraWeights->GetMCSpectraWeight(fMCLabel, 0);
        fMCweightSysUp = fMCSpectraWeights->GetMCSpectraWeight(fMCLabel, 1);
        fMCweightSysDown = fMCSpectraWeights->GetMCSpectraWeight(fMCLabel, -1);
    }
    if(fMCSpectraWeights && 1==fMCPrimSec && !fMCPileUpTrack && fMCParticle->Particle()){ // only for secondaries from decay
        fMCweight = fMCSpectraWeights->GetWeightForSecondaryParticle(fMCLabel);
        fMCweightSysUp = fMCSpectraWeights->GetWeightForSecondaryParticle(fMCLabel, 1);
        fMCweightSysDown = fMCSpectraWeights->GetWeightForSecondaryParticle(fMCLabel, -1);
    }

    for (int i = 0; i < 23; ++i) {
        if (fAcceptTrack[i]) {
            fHistEffContNCluster.FillWeight(fMCweight,fMultPercentileV0M, fPt, fTPCSignalN, i,
                                      fMCProdcutionType);
            fHistEffContZ.FillWeight(fMCweight,fMultPercentileV0M, fPt, fZInner, i,
                               fMCProdcutionType, fEta);
            fHistEffContEta.FillWeight(fMCweight,fMultPercentileV0M, fPt, fEta, i,
                                       fMCProdcutionType, fChargeSign, 0);
            fHistEffContEta.FillWeight(fMCweightSysUp,fMultPercentileV0M, fPt, fEta, i,
                                       fMCProdcutionType, fChargeSign, 1);
            fHistEffContEta.FillWeight(fMCweightSysDown,fMultPercentileV0M, fPt, fEta, i,
                                       fMCProdcutionType, fChargeSign, -1);
            fHistEffContPhi.FillWeight(fMCweight,fMultPercentileV0M, fPt, fPhi, i,
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
    if (TMath::Abs(fMCEta) > 1.0)
        return;
    if(fMCPileUpTrack){ // reject pileup tracks
        return;
    }

    if (fMCParticleType == AlidNdPtTools::kOther) {
        Log("GenPrim.PDG.", fMCPDGCode);
    }
    if (TMath::Abs(fMCQ > 1)) {
        Log("GenPrim.Q>1.PDG.", fMCPDGCode);
    }

    double fMCweight = 1.0;
    double fMCweightSysUp = 1.0;
    double fMCweightSysDown = 1.0;
    if(fMCSpectraWeights && 0==fMCPrimSec && !fMCPileUpTrack && fMCParticle->Particle()){ // only for primary particles
        fMCweight = fMCSpectraWeights->GetMCSpectraWeight(fMCLabel, 0);
        fMCweightSysUp = fMCSpectraWeights->GetMCSpectraWeight(fMCLabel, 1);
        fMCweightSysDown = fMCSpectraWeights->GetMCSpectraWeight(fMCLabel, -1);
    }
    if(fMCSpectraWeights && 1==fMCPrimSec && !fMCPileUpTrack && fMCParticle->Particle()){ // only for secondaries from decay
        fMCweight = fMCSpectraWeights->GetWeightForSecondaryParticle(fMCLabel);
        fMCweightSysUp = fMCSpectraWeights->GetWeightForSecondaryParticle(fMCLabel, 1);
        fMCweightSysDown = fMCSpectraWeights->GetWeightForSecondaryParticle(fMCLabel, -1);
    }

    fHistEffContNCluster.FillWeight(fMCweight,fMultPercentileV0M, fMCPt, fTPCSignalN, -1., 3.);
    fHistEffContZ.FillWeight(fMCweight,fMultPercentileV0M, fMCPt, fZInner, -1., 3., fMCEta);
    fHistEffContEta.FillWeight(fMCweight,fMultPercentileV0M, fMCPt, fMCEta, -1., 3., fChargeSign, 0);
    fHistEffContEta.FillWeight(fMCweightSysUp,fMultPercentileV0M, fMCPt, fMCEta, -1., 3., fChargeSign, 1);
    fHistEffContEta.FillWeight(fMCweightSysDown,fMultPercentileV0M, fMCPt, fMCEta, -1., 3., fChargeSign, -1);
    fHistEffContPhi.FillWeight(fMCweight,fMultPercentileV0M, fMCPt, fMCPhi, -1., 3.);
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
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaulteta10"));
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
        7, AlidNdPtTools::CreateESDtrackCuts("defaulteta10"));
    for(int i=1; i<16; ++i){
        task->SetESDtrackCuts(
                              i+7,AlidNdPtTools::CreateESDtrackCuts("defaultEta10", i+100));
    }

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
