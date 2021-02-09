#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskBaseWeights.h"
#include "AliAnalysisTaskMKBase.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliGenEventHeader.h"
#include "AliHeader.h"
#include "AliMCEvent.h"
#include "AliMCSpectraWeights.h"
#include "AlidNdPtTools.h"
#include "TChain.h"
#include "TGeoGlobalMagField.h"
#include "TH1F.h"
#include "TList.h"
#include "TRandom3.h"
#include <iostream>
#include <vector>

#ifdef __AliAnalysisTaskBaseWeights_DebugPCC__
#define DebugPCC(x) std::cout << x
#else
#define DebugPCC(x)
#endif

class AliAnalysisTaskBaseWeights;

namespace {
using namespace Hist;
using namespace std;
} // namespace

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskBaseWeights)
/// \endcond
//_____________________________________________________________________________

AliAnalysisTaskBaseWeights::AliAnalysisTaskBaseWeights()
: AliAnalysisTaskMKBase(), fCollSystem(CollisionSystem::pp),
fUseRandomSeed(kFALSE), fRand(0), fMCSpectraWeights(0), fMCweight(1),
fMCweightRandom(1), fMCweightSys(1), fMCweightSysRandom(1), fNch(0),
fNchWeighted(0), fNchWeightedRandom(0), fNchWeightedSys(0),
fNchWeightedSysRandom(0), fNacc(0), fNaccWeighted(0),
fNaccWeightedRandom(0), fNaccWeightedSys(0),
fNaccWeightedSysRandom(0), fHistEffCont{}, fHistMultCorrelation{} {
    // default contructor
}

//_____________________________________________________________________________

AliAnalysisTaskBaseWeights::AliAnalysisTaskBaseWeights(const char* name)
: AliAnalysisTaskMKBase(name), fCollSystem(CollisionSystem::pp),
fUseRandomSeed(kFALSE), fRand(0), fMCSpectraWeights(0), fMCweight(1),
fMCweightRandom(1), fMCweightSys(1), fMCweightSysRandom(1), fNch(0),
fNchWeighted(0), fNchWeightedRandom(0), fNchWeightedSys(0),
fNchWeightedSysRandom(0), fNacc(0), fNaccWeighted(0),
fNaccWeightedRandom(0), fNaccWeightedSys(0),
fNaccWeightedSysRandom(0), fHistEffCont{}, fHistMultCorrelation{} {
    // constructor
}
//_____________________________________________________________________________

AliAnalysisTaskBaseWeights::~AliAnalysisTaskBaseWeights() {
    if (fRand) {
        delete fRand;
        fRand = 0;
    }
}

//_____________________________________________________________________________

void AliAnalysisTaskBaseWeights::AddOutput() {
    // add also all the default output from the base class
    //    AliAnalysisTaskMKBase::BaseAddOutput();

    std::vector<double> centBins = {0.,  10., 20., 30., 40.,
        50., 60., 70., 80., 90.};
    std::vector<double> ptBins = {
        0.0,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45, 0.5,  0.55,  0.6,
        0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95, 1.0,  1.1,  1.2,   1.4,
        1.6,  1.8,  2.0,  2.2,  2.4,  2.6,  2.8,  3.0,  3.2,  3.6,   4.0,
        5.0,  6.0,  8.0,  10.0, 13.0, 20.0, 30.0, 50.0, 80.0, 100.0, 200.0};

    std::vector<double> multBins;
    if (fCollSystem == AliAnalysisTaskBaseWeights::CollisionSystem::pp) {
        multBins.reserve(51);
        for (int i = 0; i < 51; ++i) {
            multBins.push_back(i);
        }
    } else if (fCollSystem ==
               AliAnalysisTaskBaseWeights::CollisionSystem::pPb) {
        multBins.reserve(301);
        for (int i = 0; i < 301; ++i) {
            multBins.push_back(i);
        }
    } else if (fCollSystem ==
               AliAnalysisTaskBaseWeights::CollisionSystem::PbPb) {
        multBins.reserve(201);
        for (int i = 0; i < 201; ++i) {
            multBins.push_back(i * 25);
        }
    } else if (fCollSystem ==
               AliAnalysisTaskBaseWeights::CollisionSystem::XeXe) {
        multBins.reserve(201);
        for (int i = 0; i < 201; ++i) {
            multBins.push_back(i * 25);
        }
    } else {
        multBins.reserve(51);
        for (int i = 0; i < 51; ++i) {
            multBins.push_back(i);
        }
    }

    Axis centAxis{"cent", "centrality", centBins};
    Axis multAxisNch{"Nch", "#it{N}_{ch}", multBins};
    Axis multAxisNacc{"Nacc", "#it{N}_{acc}", multBins};
    Axis ptAxis{"pt", "#it{p}_{T} (GeV/c)", ptBins};
    Axis pidAxis {
        "pid", "pid", {-0.5, 9.5}, 10}; // 0=e, 1=mu, 2=pi, 3=K, 4=p, 6=sigmaP,
                                        // 7=sigmaM, 8=xi, 9=omega, 5=other
    Axis mcInfoAxis{"mcInfo",
        "mcInfo",
        {-0.3, 3.5},
        4}; // 0=prim, 1=decay 2=material, 3=genprim
    Axis mcWeightAxis{"weight",
        "weight",
        {-0.5, 4.5},
        5}; // 0=none, 1=weighted 2=weightedRandom, 3=weightSys,
            // 4=weightSysRandom

    // Hists
    double requiredMemory = 0.;
    fHistEffCont.AddAxis(centAxis);
    fHistEffCont.AddAxis(multAxisNch);
    fHistEffCont.AddAxis(ptAxis);
    fHistEffCont.AddAxis(pidAxis);
    fHistEffCont.AddAxis(mcInfoAxis);
    fHistEffCont.AddAxis(mcWeightAxis);
    fOutputList->Add(fHistEffCont.GenerateHist("fHistEffCont"));
    requiredMemory += fHistEffCont.GetSize();

    fHistMultCorrelation.AddAxis(centAxis);
    fHistMultCorrelation.AddAxis(multAxisNch);
    fHistMultCorrelation.AddAxis(multAxisNacc);
    fOutputList->Add(fHistMultCorrelation.GenerateHist("fHistMultCorrelation"));
    requiredMemory += fHistEffCont.GetSize();

    AliError(Form("Estimated memory usage of histograms: %.0f Bytes (%f MiB)",
                  requiredMemory, requiredMemory / 1048576));
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskBaseWeights::IsEventSelected() {
    return fIsAcceptedAliEventCuts;
}

//_____________________________________________________________________________

UInt_t AliAnalysisTaskBaseWeights::GetSeed() {
    if (fUseRandomSeed) {
        return 0;
    }

    UInt_t seed = fEventNumberInFile;
    seed <<= 7;
    seed += fRunNumber;
    seed <<= 7;
    seed += fMCLabel;
    seed <<= 7;
    seed += fTimeStamp;
//    DebugPCC("Seed:" << seed << "\n");
    return seed;
}

//_____________________________________________________________________________

void AliAnalysisTaskBaseWeights::AnaEventMC() {
    fNch = 0;
    fNchWeighted = 0;
    fNchWeightedRandom = 0;
    fNchWeightedSys = 0;
    fNchWeightedSysRandom = 0;
    fNacc = 0;
    fNaccWeighted = 0;
    fNaccWeightedRandom = 0;
    fNaccWeightedSys = 0;
    fNaccWeightedSysRandom = 0;

    fMCSpectraWeights = static_cast<AliMCSpectraWeightsHandler*>(
                                                                 fEvent->FindListObject("fMCSpectraWeights"))
    ->fMCSpectraWeight;

    if (fMCSpectraWeights) {
        DebugPCC("found fMCSpectraWeights in this event\n");
        DebugPCC("Status: " << fMCSpectraWeights->GetTaskStatus() << "\n");
    } else {
        DebugPCC("could not find fMCSpectraWeights in this event\n");
    }

    LoopOverAllTracks();
    LoopOverAllParticles();

    fHistMultCorrelation.Fill(static_cast<Double_t>(fNch), static_cast<Double_t>(fNacc), static_cast<Double_t>(0));
    fHistMultCorrelation.Fill(static_cast<Double_t>(fNchWeighted), static_cast<Double_t>(fNaccWeighted), static_cast<Double_t>(1));
    fHistMultCorrelation.Fill(static_cast<Double_t>(fNchWeightedRandom), static_cast<Double_t>(fNaccWeightedRandom), static_cast<Double_t>(2));
    fHistMultCorrelation.Fill(static_cast<Double_t>(fNchWeightedSys), static_cast<Double_t>(fNaccWeightedSys), static_cast<Double_t>(3));
    fHistMultCorrelation.Fill(static_cast<Double_t>(fNchWeightedSysRandom), static_cast<Double_t>(fNaccWeightedSysRandom), static_cast<Double_t>(4));

    DebugPCC("\tmultiplicities:\t mult\t weighted\t weightedRandom\t weightSys\t "
             "weightSysRandom\n");
    DebugPCC("Nch \t "
             << fNch << "\t "
             << fNchWeighted << "\t "
             << fNchWeightedRandom << "\t "
             << fNchWeightedSys << "\t "
             << fNchWeightedSysRandom << "\n\n");
    DebugPCC("Nacc \t "
             << fNacc << "\t "
             << fNaccWeighted << "\t "
             << fNaccWeightedRandom << "\t "
             << fNaccWeightedSys << "\t "
             << fNaccWeightedSysRandom << "\n");

    DebugPCC("---------- event end ---------\n");
    DebugPCC("\n\n\n");
}

//_____________________________________________________________________________

double AliAnalysisTaskBaseWeights::GetRandomRoundDouble(double val) {
    double result = 0;
    while (val >= 1) {
        ++result;
        --val;
    }
    if (val > 0) {
        if (!fRand) {
            fRand = new TRandom3();
        }
        fRand->SetSeed(GetSeed());
        if (fRand->Rndm() < val) {
            ++result;
        }
    }
    return result;
}

//_____________________________________________________________________________

void AliAnalysisTaskBaseWeights::AnaTrackMC(Int_t flag) {
    if (!fAcceptTrackM)
        return;

    if (fMCParticleType == AlidNdPtTools::kOther) {
        Log("RecTrack.PDG.", fMCPDGCode);
    }
    if (TMath::Abs(fMCQ > 1)) {
        Log("RecTrack.Q>1.PDG.", fMCPDGCode);
    }

    fHistEffCont.FillWeight(
                            static_cast<Double_t>(1), static_cast<Double_t>(fMultPercentileV0M),
                            static_cast<Double_t>(fNTracksAcc), static_cast<Double_t>(fMCPt),
                            static_cast<Double_t>(fMCParticleType),
                            static_cast<Double_t>(fMCProdcutionType), static_cast<Double_t>(0));
    fHistEffCont.FillWeight(
                            fMCweight, fMultPercentileV0M, static_cast<Double_t>(fNTracksAcc),
                            static_cast<Double_t>(fMCPt), static_cast<Double_t>(fMCParticleType),
                            static_cast<Double_t>(fMCProdcutionType), static_cast<Double_t>(1));
    fHistEffCont.FillWeight(
                            fMCweightRandom, static_cast<Double_t>(fMultPercentileV0M),
                            static_cast<Double_t>(fNTracksAcc), fMCPt,
                            static_cast<Double_t>(fMCParticleType),
                            static_cast<Double_t>(fMCProdcutionType), static_cast<Double_t>(2));
    fHistEffCont.FillWeight(
                            fMCweightSys, static_cast<Double_t>(fMultPercentileV0M),
                            static_cast<Double_t>(fNTracksAcc), static_cast<Double_t>(fMCPt),
                            static_cast<Double_t>(fMCParticleType),
                            static_cast<Double_t>(fMCProdcutionType), static_cast<Double_t>(3));
    fHistEffCont.FillWeight(
                            fMCweightSysRandom, static_cast<Double_t>(fMultPercentileV0M),
                            static_cast<Double_t>(fNTracksAcc), static_cast<Double_t>(fMCPt),
                            static_cast<Double_t>(fMCParticleType),
                            static_cast<Double_t>(fMCProdcutionType), static_cast<Double_t>(4));

    if (fPt > 0.15 && fPt < 50.0 && TMath::Abs(fEta) < 0.8) {
        ++fNacc;
        fNaccWeighted += fMCweight;
        fNaccWeightedRandom += fMCweightRandom;
        fNaccWeightedSys += fMCweightSys;
        fNaccWeightedSysRandom += fMCweightSysRandom;
    }
}

//_____________________________________________________________________________

void AliAnalysisTaskBaseWeights::AnaParticleMC(Int_t flag) {
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

    fHistEffCont.FillWeight(
                             static_cast<Double_t>(1), fMultPercentileV0M, static_cast<Double_t>(fNTracksAcc), fMCPt,
                             static_cast<Double_t>(fMCParticleType), static_cast<Double_t>(3), static_cast<Double_t>(0));
    fHistEffCont.FillWeight(
                             static_cast<Double_t>(fMCweight), fMultPercentileV0M, static_cast<Double_t>(fNTracksAcc), fMCPt,
                             static_cast<Double_t>(fMCParticleType), static_cast<Double_t>(3), static_cast<Double_t>(1));
    fHistEffCont.FillWeight(
                             static_cast<Double_t>(fMCweightRandom), static_cast<Double_t>(fMultPercentileV0M), static_cast<Double_t>(fNTracksAcc),
                             static_cast<Double_t>(fMCPt), static_cast<Double_t>(fMCParticleType), static_cast<Double_t>(3), static_cast<Double_t>(2));
    fHistEffCont.FillWeight(
                             static_cast<Double_t>(fMCweightSys), static_cast<Double_t>(fMultPercentileV0M), static_cast<Double_t>(fNTracksAcc),
                             static_cast<Double_t>(fMCPt), static_cast<Double_t>(fMCParticleType), static_cast<Double_t>(3), static_cast<Double_t>(3));
    fHistEffCont.FillWeight(static_cast<Double_t>(fMCweightSysRandom), static_cast<Double_t>(fMultPercentileV0M),
                             static_cast<Double_t>(fNTracksAcc), static_cast<Double_t>(fMCPt),
                             static_cast<Double_t>(fMCParticleType),
                             static_cast<Double_t>(3), static_cast<Double_t>(4));

    if (fPt > 0.15 && fPt < 50.0 && TMath::Abs(fEta) < 0.8) {
        ++fNch;
        fNchWeighted += fMCweight;
        fNchWeightedRandom += fMCweightRandom;
        fNchWeightedSys += fMCweightSys;
        fNchWeightedSysRandom += fMCweightSysRandom;
    }
}

//_____________________________________________________________________________

void AliAnalysisTaskBaseWeights::LoopOverAllTracks(Int_t flag) {
    // on data or if fUseMCWeights is not set we do not modify anything
    if (!fIsMC) {
        AliAnalysisTaskMKBase::LoopOverAllTracks(flag);
    }

    fNTracksESD = fESD->GetNumberOfTracks();
    for (Int_t i = 0; i < fNTracksESD; i++) {
        fESDTrack = dynamic_cast<AliESDtrack*>(fESD->GetTrack(i));
        if (!fESDTrack) {
            Err("noESDtrack");
            continue;
        }
        InitTrack();

        // get the scaling factor
        fMCweight = fMCSpectraWeights->GetMCSpectraWeightNominal(
                                                                 fMCParticle->Particle());
        fMCweightSys = fMCSpectraWeights->GetMCSpectraWeightSystematics(
                                                                        fMCParticle->Particle());
        fMCweightRandom = GetRandomRoundDouble(fMCweight);
        fMCweightSysRandom = GetRandomRoundDouble(fMCweightSys);
        BaseAnaTrack(flag);
    }
}

//_____________________________________________________________________________

void AliAnalysisTaskBaseWeights::LoopOverAllParticles(Int_t flag) {
    // this method should not be called on data
    if (!fIsMC)
        return;

    fMCnTracks = fMC->GetNumberOfTracks();
    for (Int_t i = 0; i < fMCnTracks; i++) {
        fMCParticle = dynamic_cast<AliMCParticle*>(fMC->GetTrack(i));
        if (!fMCParticle) {
            Err("noMCParticle");
            continue;
        }
        fMCLabel = i;
        InitMCParticle();
        // get the scaling factor
        fMCweight = fMCSpectraWeights->GetMCSpectraWeightNominal(
                                                                 fMCParticle->Particle());
        fMCweightSys = fMCSpectraWeights->GetMCSpectraWeightSystematics(
                                                                        fMCParticle->Particle());
        fMCweightRandom = GetRandomRoundDouble(fMCweight);
        fMCweightSysRandom = GetRandomRoundDouble(fMCweightSys);
        BaseAnaParticleMC(flag);
    }
}

//_____________________________________________________________________________

Double_t AliAnalysisTaskBaseWeights::MCScalingFactor() {
    // determine the MC scaling factor for the current particle

    // in case mcspectraweights are there we use them for primary particles
    if (fMCSpectraWeights && fMCisPrim) {
        return fMCSpectraWeights->GetMCSpectraWeight(fMCParticle->Particle(),
                                                     fMC);
    }
    if (!fMCisPrim) {
        // TODO: write secondary scaling interface here
    }

    return 1.0;
}

//_____________________________________________________________________________

AliAnalysisTaskBaseWeights* AliAnalysisTaskBaseWeights::AddTaskBaseWeights(
                                                                           const char* name, const char* outfile, CollisionSystem collisionSystem) {
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskBaseWeights", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the
    // analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskBaseWeights",
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
    AliAnalysisTaskBaseWeights* task = new AliAnalysisTaskBaseWeights(name);
    if (!task) {
        return 0;
    }

    // configure the task
    //===========================================================================
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
    task->SetESDtrackCuts(0, AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
    task->SetNeedEventMult(kTRUE);
    task->fCollSystem = collisionSystem;

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
