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

class AliAnalysisTaskBaseWeights;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskBaseWeights)
    /// \endcond
    //_____________________________________________________________________________

    AliAnalysisTaskBaseWeights::AliAnalysisTaskBaseWeights()
    : AliAnalysisTaskMKBase(), fUseMCWeights(kTRUE), fUseRandomSeed(kFALSE),
      fRand(0), fMCSpectraWeights(0), fMCweight(1),fMCweightRandom(1), fNch(0), fNchWeighted(0), fNchWeightedRandom(0), fNacc(0), fNaccWeighted(0), fNaccWeightedRandom(0), fHistEffCont(0),
      fHistMultCorrelation(0) {
    // default contructor
}

//_____________________________________________________________________________

AliAnalysisTaskBaseWeights::AliAnalysisTaskBaseWeights(const char* name)
    : AliAnalysisTaskMKBase(name), fUseMCWeights(kTRUE), fUseRandomSeed(kFALSE),
      fRand(0), fMCSpectraWeights(0), fMCweight(1),fMCweightRandom(1), fNch(0), fNchWeighted(0), fNchWeightedRandom(0), fNacc(0), fNaccWeighted(0), fNaccWeightedRandom(0), fHistEffCont(0),
      fHistMultCorrelation(0) {
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
    AliAnalysisTaskMKBase::BaseAddOutput();

    AddAxis("cent");
    AddAxis("nAcc", "mult6kcoarse");
    AddAxis("MCpT", "pt");
    AddAxis("MCQ", 3, -1.5, 1.5);
    AddAxis("MCpid", 10, -0.5, 9.5); // 0=e, 1=mu, 2=pi, 3=K, 4=p, 6=sigmaP,
                                     // 7=sigmaM, 8=xi, 9=omega, 5=other
    AddAxis("MCinfo", 4, -0.5, 3.5); // 0=prim, 1=decay 2=material, 3=genprim
    AddAxis("MCWeighted", 3, -0.5, 2.5); // 0=none, 1=weighted 2=weightedRandom
    fHistEffCont = CreateHist("fHistEffCont");
    fOutputList->Add(fHistEffCont);

    AddAxis("nAcc", "mult6kfine");
    AddAxis("nCh", "mult6kfine");
    AddAxis("MCWeighted", 3, -0.5, 2.5); // 0=none, 1=weighted 2=weightedRandom
    fHistMultCorrelation = CreateHist("fHistMultCorrelation");
    fOutputList->Add(fHistMultCorrelation);

    // if there are weights and they are used add them to the ouput
    if (fMCSpectraWeights && fUseMCWeights) {
        fOutputList->Add((TObject*)fMCSpectraWeights->GetHistMCGenPrimTrackParticles());
        fOutputList->Add((TObject*)fMCSpectraWeights->GetHistDataFraction());
        fOutputList->Add((TObject*)fMCSpectraWeights->GetHistMCFraction());
        fOutputList->Add((TObject*)fMCSpectraWeights->GetHistMCWeights());
    }
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskBaseWeights::IsEventSelected()
{
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
    return seed;
}

//_____________________________________________________________________________

void AliAnalysisTaskBaseWeights::AnaEvent()
{
   fNch=0;  fNchWeighted=0;     fNchWeightedRandom=0;
   fNacc=0; fNaccWeighted=0;    fNaccWeightedRandom=0;
    
   LoopOverAllTracks();
   if (fIsMC) LoopOverAllParticles();
    
    FillHistWeighted(fHistMultCorrelation, {fNch, fNacc, 0}, 1);
    FillHistWeighted(fHistMultCorrelation, {fNchWeighted, fNaccWeighted, 1}, 1);
    FillHistWeighted(fHistMultCorrelation, {fNchWeightedRandom, fNaccWeightedRandom, 2}, 1);
    
    std::cout << "multiplicities:\t N_ch\t weighted\t weightedRandom\n";
    std::cout << "\t " << fNch << "\t " << fNchWeighted << "\t " << fNchWeightedRandom << "\n";
    std::cout << "\t " << fNacc << "\t " << fNaccWeighted << "\t " << fNaccWeightedRandom << "\n";
}

//_____________________________________________________________________________

double AliAnalysisTaskBaseWeights::GetRandomRoundDouble(double val){
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

void AliAnalysisTaskBaseWeights::AnaTrackMC(Int_t flag)
{
    if (!fAcceptTrackM) return;
     
    if (fMCParticleType==AlidNdPtTools::kOther) { Log("RecTrack.PDG.",fMCPDGCode); }
    if (TMath::Abs(fMCQ > 1)) { Log("RecTrack.Q>1.PDG.",fMCPDGCode); }
    
    
    FillHistWeighted(fHistEffCont, {fMultPercentileV0M, static_cast<double>(fNTracksAcc), fMCPt, static_cast<double>(fMCChargeSign), static_cast<double>(fMCParticleType), static_cast<double>(fMCProdcutionType), 0}, 1);
    FillHistWeighted(fHistEffCont, {fMultPercentileV0M, static_cast<double>(fNTracksAcc), fMCPt, static_cast<double>(fMCChargeSign), static_cast<double>(fMCParticleType), static_cast<double>(fMCProdcutionType), 1}, fMCweight);
    FillHistWeighted(fHistEffCont, {fMultPercentileV0M, static_cast<double>(fNTracksAcc), fMCPt, static_cast<double>(fMCChargeSign), static_cast<double>(fMCParticleType), static_cast<double>(fMCProdcutionType), 2}, fMCweightRandom);
    
    if (fPt > 0.15 && fPt < 50.0 && TMath::Abs(fEta) < 0.8 ) {
        ++fNacc;
        fNaccWeighted=fNaccWeighted+fMCweight;
        fNaccWeightedRandom=fNaccWeightedRandom+fMCweightRandom;
    }
}

//_____________________________________________________________________________

void AliAnalysisTaskBaseWeights::AnaParticleMC(Int_t flag)
{
    if (!fMCisPrim) return;
    if (!fMCIsCharged) return;
    if (TMath::Abs(fMCEta) > 0.8) return;
    
    if (fMCParticleType==AlidNdPtTools::kOther) { Log("GenPrim.PDG.",fMCPDGCode); }
    if (TMath::Abs(fMCQ > 1)) { Log("GenPrim.Q>1.PDG.",fMCPDGCode); }
    
    
    FillHistWeighted(fHistEffCont, {fMultPercentileV0M, static_cast<double>(fNTracksAcc), fMCPt, static_cast<double>(fMCChargeSign), static_cast<double>(fMCParticleType), 3.0, 0.0}, 1);
    FillHistWeighted(fHistEffCont, {fMultPercentileV0M, static_cast<double>(fNTracksAcc), fMCPt, static_cast<double>(fMCChargeSign), static_cast<double>(fMCParticleType), 3.0, 1.0}, fMCweight);
    FillHistWeighted(fHistEffCont, {fMultPercentileV0M, static_cast<double>(fNTracksAcc), fMCPt, static_cast<double>(fMCChargeSign), static_cast<double>(fMCParticleType), 3.0, 2.0}, fMCweightRandom);
    
    if (fPt > 0.15 && fPt < 50.0 && TMath::Abs(fEta) < 0.8 ) {
        ++fNch;
        fNchWeighted=fNchWeighted+fMCweight;
        fNchWeightedRandom=fNchWeightedRandom+fMCweightRandom;
    }
}

//_____________________________________________________________________________

void AliAnalysisTaskBaseWeights::LoopOverAllTracks(Int_t flag) {
    // on data or if fUseMCWeights is not set we do not modify anything
    if (!fIsMC) {
        AliAnalysisTaskMKBase::LoopOverAllTracks(flag);
    }
    if (!fUseMCWeights) {
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
        fMCweight = MCScalingFactor();
        fMCweightRandom = GetRandomRoundDouble(fMCweight);
        BaseAnaTrack(flag);
    }
}

//_____________________________________________________________________________

void AliAnalysisTaskBaseWeights::LoopOverAllParticles(Int_t flag) {
    // this method should not be called on data
    if (!fIsMC)
        return;
    // only if fUseMCWeights is set we do re weighting
    if (!fUseMCWeights) {
        AliAnalysisTaskMKBase::LoopOverAllParticles(flag);
    }

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
        fMCweight = MCScalingFactor();
        fMCweightRandom = GetRandomRoundDouble(fMCweight);
        BaseAnaParticleMC(flag);
    }
}

//_____________________________________________________________________________

Double_t AliAnalysisTaskBaseWeights::MCScalingFactor() {
    // determine the MC scaling factor for the current particle

    // in case the weights are not activated, they are always unity
    if (!fUseMCWeights) {
        return 1.0;
    }

    // in case mcspectraweights are there we use them for primary particles
    if (fMCSpectraWeights && fMCisPrim) {
        return fMCSpectraWeights->GetMCSpectraWeight(fMCParticle->Particle(),
                                                     fMC);
    }
    if(!fMCisPrim){
    // TODO: write secondary scaling interface here
    }
    
    return 1.0;
}

//_____________________________________________________________________________

void AliAnalysisTaskBaseWeights::FillDefaultHistograms(Int_t step) {
    // fill the spectra weights if needed
    // but do it before any event selection
    // because this might depent on the derived task

    // protection
    if (fMCSpectraWeights && step == 1) {
        if (fMCSpectraWeights->GetTaskStatus() <
            AliMCSpectraWeights::TaskState::kMCSpectraObtained) {
            // for now I pass the V0M multiplicity percentile
            fMCSpectraWeights->FillMCSpectra(fMC);
            //            cout<<"fMCSpectraWeights->FillMCSpectra(fMC);"<<endl;
            //            //DEBUG
        }
    }

    // also fill the default histograms
    AliAnalysisTaskMKBase::FillDefaultHistograms(step);
}

//_____________________________________________________________________________

AliAnalysisTaskBaseWeights* AliAnalysisTaskBaseWeights::AddTaskBaseWeights(
    const char* name, const char* outfile, const char* collisionSystem,
    Int_t sysFlag, const char* prevTrainOutputPath) {
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
    // configure the use of AliMCSpectraWeights
    //===========================================================================
    // collisionSystem is String "pp", "pPb", "XeXe", "PbPb"
    if (collisionSystem) {
        AliMCSpectraWeights* weights =
            new AliMCSpectraWeights(collisionSystem, "fMCSpectraWeights",
                                    (AliMCSpectraWeights::SysFlag)sysFlag);
        // root file with fHistMCGenPrimTrackParticle (THnF)
        // or used SetSavedListName
        // prevTrainOutputPath is string, for lego train activate the check box
        // addtask to get alien connection for addtask file is only needed for
        // addtask
        //
        // file with fractions from data: expert input from Patrick:
        // /alice/cern.ch/user/p/phuhn/AllPublishedFractions.root
        // (path fixed in AliMCSpectraWeights code)
        //  used void SetDataFractionsFile(const char* file) to set local path

        if (prevTrainOutputPath) {
            weights->SetMCSpectraFile(prevTrainOutputPath);
        } // path to previous train output
        weights->Init();
        task->fMCSpectraWeights = weights;
    }

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
