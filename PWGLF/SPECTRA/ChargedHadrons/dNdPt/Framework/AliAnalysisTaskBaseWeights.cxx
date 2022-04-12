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

using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskBaseWeights)
/// \endcond
//_____________________________________________________________________________

AliAnalysisTaskBaseWeights::AliAnalysisTaskBaseWeights()
: AliAnalysisTaskMKBase(), fCollSystem(CollisionSystem::pp),
fUseRandomSeed(kFALSE), fRand(0), fMCSpectraWeights(0), fMCweight(1),
fMCweightRandom(1), fMCweightSys(1), fMCweightSysDown(1), fMCweightSysRandom(1), fNch(0),
fNchWeighted(0), fNchWeightedRandom(0), fNchWeightedSys(0),
fNchWeightedSysRandom(0), fNacc(0), fNaccWeighted(0),
fNaccWeightedRandom(0), fNaccWeightedSys(0),
fNaccWeightedSysRandom(0), fHistEffContNominal{}, fHistEffContWeighted{}, fHistEffContWeightedRandom{}, fHistEffContWeightedSys{}, fHistEffContWeightedSysDown{}, fHistEffContWeightedSysRandom{}, fHistMultCorrelationNominal{}, fHistMultCorrelationWeighted{}, fHistMultCorrelationWeightedRandom{}, fHistMultCorrelationWeightedSys{}, fHistMultCorrelationWeightedSysRandom{}, fHistPionRec{}, fHistPionGen{} {
    // default contructor
}

//_____________________________________________________________________________

AliAnalysisTaskBaseWeights::AliAnalysisTaskBaseWeights(const char* name)
: AliAnalysisTaskMKBase(name), fCollSystem(CollisionSystem::pp),
fUseRandomSeed(kFALSE), fRand(0), fMCSpectraWeights(0), fMCweight(1),
fMCweightRandom(1), fMCweightSys(1),fMCweightSysDown(1), fMCweightSysRandom(1), fNch(0),
fNchWeighted(0), fNchWeightedRandom(0), fNchWeightedSys(0),
fNchWeightedSysRandom(0), fNacc(0), fNaccWeighted(0),
fNaccWeightedRandom(0), fNaccWeightedSys(0),
fNaccWeightedSysRandom(0), fHistEffContNominal{}, fHistEffContWeighted{}, fHistEffContWeightedRandom{}, fHistEffContWeightedSys{}, fHistEffContWeightedSysDown{}, fHistEffContWeightedSysRandom{}, fHistMultCorrelationNominal{}, fHistMultCorrelationWeighted{}, fHistMultCorrelationWeightedRandom{}, fHistMultCorrelationWeightedSys{}, fHistMultCorrelationWeightedSysRandom{}, fHistPionRec{}, fHistPionGen{} {
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
    std::vector<double> ptBins = {0.0, 0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.5,5,5.5,6,6.5,7,8,10,13,20, 30, 50, 80, 100, 200};

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
    std::vector<double> etaBins{-1,-0.8,0.8,1};
    std::vector<double> yBins{-1,-0.5,0.5,1};

    Axis etaAxis{"#eta", "pseudorapidity", etaBins};
    Axis yAxis{"y", "rapidity", yBins};
    Axis centAxis{"cent", "centrality", centBins};
    Axis multAxisNch{"Nch", "#it{N}_{ch}", multBins};
    Axis multAxisNacc{"Nacc", "#it{N}_{acc}", multBins};
    Axis ptAxis{"pt", "#it{p}_{T} (GeV/c)", ptBins};
    Axis pidAxis {
        "pid", "pid", {-0.5, 9.5}, 10}; // 0=e, 1=mu, 2=pi, 3=K, 4=p, 6=sigmaP,
                                        // 7=sigmaM, 8=xi, 9=omega, 5=other
    Axis mcInfoAxis{"mcInfo",
        "mcInfo",
        {-0.5, 3.5},
        4}; // 0=prim, 1=decay 2=material, 3=genprim
    Axis mcWeightAxis{"weight",
        "weight",
        {-0.5, 4.5},
        5}; // 0=none, 1=weighted 2=weightedRandom, 3=weightSys,
            // 4=weightSysRandom

    // Hists
    double requiredMemory = 0.;
//    fHistEffCont.AddAxis(centAxis);
    fHistEffContNominal.AddAxis(multAxisNch);
    fHistEffContNominal.AddAxis(ptAxis);
    fHistEffContNominal.AddAxis(pidAxis);
    fHistEffContNominal.AddAxis(mcInfoAxis);
//    fHistEffCont.AddAxis(mcWeightAxis);
    fHistEffContNominal.AddAxis(etaAxis);
    fHistEffContNominal.AddAxis(yAxis);
    fOutputList->Add(fHistEffContNominal.GenerateHist("fHistEffContNominal"));
    requiredMemory += fHistEffContNominal.GetSize();

    fHistEffContWeighted.AddAxis(multAxisNch);
    fHistEffContWeighted.AddAxis(ptAxis);
    fHistEffContWeighted.AddAxis(pidAxis);
    fHistEffContWeighted.AddAxis(mcInfoAxis);
    fHistEffContWeighted.AddAxis(etaAxis);
    fHistEffContWeighted.AddAxis(yAxis);
    fOutputList->Add(fHistEffContWeighted.GenerateHist("fHistEffContWeighted"));
    requiredMemory += fHistEffContWeighted.GetSize();

    fHistEffContWeightedRandom.AddAxis(multAxisNch);
    fHistEffContWeightedRandom.AddAxis(ptAxis);
    fHistEffContWeightedRandom.AddAxis(pidAxis);
    fHistEffContWeightedRandom.AddAxis(mcInfoAxis);
    fHistEffContWeightedRandom.AddAxis(etaAxis);
    fHistEffContWeightedRandom.AddAxis(yAxis);
    fOutputList->Add(fHistEffContWeightedRandom.GenerateHist("fHistEffContWeightedRandom"));
    requiredMemory += fHistEffContWeightedRandom.GetSize();

    fHistEffContWeightedSys.AddAxis(multAxisNch);
    fHistEffContWeightedSys.AddAxis(ptAxis);
    fHistEffContWeightedSys.AddAxis(pidAxis);
    fHistEffContWeightedSys.AddAxis(mcInfoAxis);
    fHistEffContWeightedSys.AddAxis(etaAxis);
    fHistEffContWeightedSys.AddAxis(yAxis);
    fOutputList->Add(fHistEffContWeightedSys.GenerateHist("fHistEffContWeightedSys"));
    requiredMemory += fHistEffContWeightedSys.GetSize();

    fHistEffContWeightedSysDown.AddAxis(multAxisNch);
    fHistEffContWeightedSysDown.AddAxis(ptAxis);
    fHistEffContWeightedSysDown.AddAxis(pidAxis);
    fHistEffContWeightedSysDown.AddAxis(mcInfoAxis);
    fHistEffContWeightedSysDown.AddAxis(etaAxis);
    fHistEffContWeightedSysDown.AddAxis(yAxis);
    fOutputList->Add(fHistEffContWeightedSysDown.GenerateHist("fHistEffContWeightedSysDown"));
    requiredMemory += fHistEffContWeightedSysDown.GetSize();

    fHistEffContWeightedSysRandom.AddAxis(multAxisNch);
    fHistEffContWeightedSysRandom.AddAxis(ptAxis);
    fHistEffContWeightedSysRandom.AddAxis(pidAxis);
    fHistEffContWeightedSysRandom.AddAxis(mcInfoAxis);
    fHistEffContWeightedSysRandom.AddAxis(etaAxis);
    fHistEffContWeightedSysRandom.AddAxis(yAxis);
    fOutputList->Add(fHistEffContWeightedSysRandom.GenerateHist("fHistEffContWeightedSysRandom"));
    requiredMemory += fHistEffContWeightedSysRandom.GetSize();

    // correlation hists
    fHistMultCorrelationNominal.AddAxis(multAxisNch);
    fHistMultCorrelationNominal.AddAxis(multAxisNacc);
    fOutputList->Add(fHistMultCorrelationNominal.GenerateHist("fHistMultCorrelationNominal"));
    requiredMemory += fHistMultCorrelationNominal.GetSize();

    fHistMultCorrelationWeighted.AddAxis(multAxisNch);
    fHistMultCorrelationWeighted.AddAxis(multAxisNacc);
    fOutputList->Add(fHistMultCorrelationWeighted.GenerateHist("fHistMultCorrelationWeighted"));
    requiredMemory += fHistMultCorrelationWeighted.GetSize();

    fHistMultCorrelationWeightedRandom.AddAxis(multAxisNch);
    fHistMultCorrelationWeightedRandom.AddAxis(multAxisNacc);
    fOutputList->Add(fHistMultCorrelationWeightedRandom.GenerateHist("fHistMultCorrelationWeightedRandom"));
    requiredMemory += fHistMultCorrelationWeightedRandom.GetSize();

    fHistMultCorrelationWeightedSys.AddAxis(multAxisNch);
    fHistMultCorrelationWeightedSys.AddAxis(multAxisNacc);
    fOutputList->Add(fHistMultCorrelationWeightedSys.GenerateHist("fHistMultCorrelationWeightedSys"));
    requiredMemory += fHistMultCorrelationWeightedSys.GetSize();

    fHistMultCorrelationWeightedSysRandom.AddAxis(multAxisNch);
    fHistMultCorrelationWeightedSysRandom.AddAxis(multAxisNacc);
    fOutputList->Add(fHistMultCorrelationWeightedSysRandom.GenerateHist("fHistMultCorrelationWeightedSysRandom"));
    requiredMemory += fHistMultCorrelationWeightedSysRandom.GetSize();

    // sanity QA
    fHistPionRec.AddAxis(multAxisNch);
    fHistPionRec.AddAxis(ptAxis);
    fOutputList->Add(fHistPionRec.GenerateHist("fHistPionRec"));
    requiredMemory += fHistPionRec.GetSize();

    fHistPionGen.AddAxis(multAxisNch);
    fHistPionGen.AddAxis(ptAxis);
    fOutputList->Add(fHistPionGen.GenerateHist("fHistPionGen"));
    requiredMemory += fHistPionGen.GetSize();

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

    AliMCSpectraWeightsHandler* mcWeightsHandler = static_cast<AliMCSpectraWeightsHandler*>(fEvent->FindListObject("fMCSpectraWeights"));
    fMCSpectraWeights = (mcWeightsHandler) ? mcWeightsHandler->fMCSpectraWeight : nullptr;

//    if (fMCSpectraWeights) {
//        DebugPCC("found fMCSpectraWeights in this event\n");
//        DebugPCC("Status: " << fMCSpectraWeights->GetTaskStatus() << "\n");
//    } else {
//        DebugPCC("could not find fMCSpectraWeights in this event\n");
//    }

    LoopOverAllParticles();
    LoopOverAllTracks();

    fHistMultCorrelationNominal.Fill(static_cast<Double_t>(fNch), static_cast<Double_t>(fNacc));
    fHistMultCorrelationWeighted.Fill(static_cast<Double_t>(fNchWeighted), static_cast<Double_t>(fNaccWeighted));
    fHistMultCorrelationWeightedRandom.Fill(static_cast<Double_t>(fNchWeightedRandom), static_cast<Double_t>(fNaccWeightedRandom));
    fHistMultCorrelationWeightedSys.Fill(static_cast<Double_t>(fNchWeightedSys), static_cast<Double_t>(fNaccWeightedSys));
    fHistMultCorrelationWeightedSysRandom.Fill(static_cast<Double_t>(fNchWeightedSysRandom), static_cast<Double_t>(fNaccWeightedSysRandom));

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
    if(fMCPileUpTrack){
        Log("PileUpTrack");
        return;
    }

    // get the scaling factor
    if(fMCSpectraWeights && 0==fMCPrimSec && fMCParticle->Particle()){ // only for primary particles
        fMCweight = fMCSpectraWeights->GetMCSpectraWeight(fMCLabel, 0);
        fMCweightSys = fMCSpectraWeights->GetMCSpectraWeight(fMCLabel, 1);
        fMCweightSysDown = fMCSpectraWeights->GetMCSpectraWeight(fMCLabel, -1);
        fMCweightRandom = GetRandomRoundDouble(fMCweight);
        fMCweightSysRandom = GetRandomRoundDouble(fMCweightSys);
    }
    if(fMCSpectraWeights && 1==fMCPrimSec && fMCParticle->Particle()){ // only for secondaries from decay
        fMCweight = fMCSpectraWeights->GetWeightForSecondaryParticle(fMCLabel);
        fMCweightSys = fMCSpectraWeights->GetWeightForSecondaryParticle(fMCLabel, 1);
        fMCweightSysDown = fMCSpectraWeights->GetWeightForSecondaryParticle(fMCLabel, -1);
        fMCweightRandom = GetRandomRoundDouble(fMCweight);
        fMCweightSysRandom = GetRandomRoundDouble(fMCweightSys);
    }

    // Fill Histograms
    fHistEffContNominal.FillWeight(
                            static_cast<Double_t>(1),
                            static_cast<Double_t>(fMCnPrim05), static_cast<Double_t>(fMCPt),
                            static_cast<Double_t>(fMCParticleType), static_cast<Double_t>(fMCProdcutionType), fMCEta, fMCY);
    fHistEffContWeighted.FillWeight(
                            fMCweight, static_cast<Double_t>(fMCnPrim05),
                            static_cast<Double_t>(fMCPt), static_cast<Double_t>(fMCParticleType),
                            static_cast<Double_t>(fMCProdcutionType), fMCEta, fMCY);
    fHistEffContWeightedRandom.FillWeight(
                            fMCweightRandom,
                            static_cast<Double_t>(fMCnPrim05), static_cast<Double_t>(fMCPt),
                            static_cast<Double_t>(fMCParticleType),
                            static_cast<Double_t>(fMCProdcutionType), fMCEta, fMCY);
    fHistEffContWeightedSys.FillWeight(
                            fMCweightSys,
                            static_cast<Double_t>(fMCnPrim05), static_cast<Double_t>(fMCPt),
                            static_cast<Double_t>(fMCParticleType),
                            static_cast<Double_t>(fMCProdcutionType), fMCEta, fMCY);
    fHistEffContWeightedSysDown.FillWeight(
                                           fMCweightSysDown, static_cast<Double_t>(fMCnPrim05),
                                           static_cast<Double_t>(fMCPt), static_cast<Double_t>(fMCParticleType),
                                           static_cast<Double_t>(fMCProdcutionType), fMCEta, fMCY);
    fHistEffContWeightedSysRandom.FillWeight(
                            fMCweightSysRandom,
                            static_cast<Double_t>(fMCnPrim05), static_cast<Double_t>(fMCPt),
                            static_cast<Double_t>(fMCParticleType),
                            static_cast<Double_t>(fMCProdcutionType), fMCEta, fMCY);

    if(fMCParticleType==AlidNdPtTools::ParticleType::kPi){
        // fill rec
        fHistPionRec.Fill(fMCnPrim05, fMCPt);
    }

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
    if (TMath::Abs(fMCEta) > 1)
        return;

    if (fMCParticleType == AlidNdPtTools::kOther) {
        Log("GenPrim.PDG.", fMCPDGCode);
    }
    if (TMath::Abs(fMCQ > 1)) {
        Log("GenPrim.Q>1.PDG.", fMCPDGCode);
    }
    if(fMCPileUpTrack){
        Log("PileUpTrack");
        return;
    }

    // get the scaling factor
    if(fMCSpectraWeights && 0==fMCPrimSec && fMCParticle->Particle()){ // only for primary particles
        fMCweight = fMCSpectraWeights->GetMCSpectraWeight(fMCLabel, 0);
        fMCweightSys = fMCSpectraWeights->GetMCSpectraWeight(fMCLabel, 1);
        fMCweightSysDown = fMCSpectraWeights->GetMCSpectraWeight(fMCLabel, -1);
        fMCweightRandom = GetRandomRoundDouble(fMCweight);
        fMCweightSysRandom = GetRandomRoundDouble(fMCweightSys);
    }
    if(fMCSpectraWeights && 1==fMCPrimSec && fMCParticle->Particle()){ // only for secondaries from decay
        fMCweight = fMCSpectraWeights->GetWeightForSecondaryParticle(fMCLabel);
        fMCweightSys = fMCSpectraWeights->GetWeightForSecondaryParticle(fMCLabel, 1);
        fMCweightSysDown = fMCSpectraWeights->GetWeightForSecondaryParticle(fMCLabel, -1);
        fMCweightRandom = GetRandomRoundDouble(fMCweight);
        fMCweightSysRandom = GetRandomRoundDouble(fMCweightSys);
    }
    // fill histograms
    fHistEffContNominal.FillWeight(
                             static_cast<Double_t>(1), static_cast<Double_t>(fMCnPrim05), fMCPt,
                             static_cast<Double_t>(fMCParticleType), static_cast<Double_t>(3), fMCEta, fMCY);
    fHistEffContWeighted.FillWeight(
                             static_cast<Double_t>(fMCweight), static_cast<Double_t>(fMCnPrim05), fMCPt,
                             static_cast<Double_t>(fMCParticleType), static_cast<Double_t>(3), fMCEta, fMCY);
    fHistEffContWeightedRandom.FillWeight(
                             static_cast<Double_t>(fMCweightRandom), static_cast<Double_t>(fMCnPrim05),
                             static_cast<Double_t>(fMCPt), static_cast<Double_t>(fMCParticleType), static_cast<Double_t>(3), fMCEta, fMCY);
    fHistEffContWeightedSys.FillWeight(
                             static_cast<Double_t>(fMCweightSys), static_cast<Double_t>(fMCnPrim05),
                             static_cast<Double_t>(fMCPt), static_cast<Double_t>(fMCParticleType), static_cast<Double_t>(3), fMCEta, fMCY);
    fHistEffContWeightedSysDown.FillWeight(
                                           static_cast<Double_t>(fMCweightSysDown), static_cast<Double_t>(fMCnPrim05),
                                           static_cast<Double_t>(fMCPt), static_cast<Double_t>(fMCParticleType), static_cast<Double_t>(3), fMCEta, fMCY);
    fHistEffContWeightedSysRandom.FillWeight(static_cast<Double_t>(fMCweightSysRandom),
                             static_cast<Double_t>(fMCnPrim05), static_cast<Double_t>(fMCPt),
                             static_cast<Double_t>(fMCParticleType),
                             static_cast<Double_t>(3), fMCEta, fMCY);

    if(fMCParticleType==AlidNdPtTools::ParticleType::kPi){
        // fill rec
        fHistPionGen.Fill(fMCnPrim05, fMCPt);
    }

    if (fPt > 0.15 && fPt < 50.0 && TMath::Abs(fEta) < 0.8) {
        ++fNch;
        fNchWeighted += fMCweight;
        fNchWeightedRandom += fMCweightRandom;
        fNchWeightedSys += fMCweightSys;
        fNchWeightedSysRandom += fMCweightSysRandom;
    }
}

//_____________________________________________________________________________

Double_t AliAnalysisTaskBaseWeights::MCScalingFactor() {
    // determine the MC scaling factor for the current particle

    // in case mcspectraweights are there we use them for primary particles
    if (fMCSpectraWeights && fMCisPrim) {
        return fMCSpectraWeights->GetMCSpectraWeight(fMCLabel,
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
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta10"));
    task->SetESDtrackCuts(0, AlidNdPtTools::CreateESDtrackCuts("defaultEta10"));
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
