#include <iostream>
#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
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
#include "AliMCSpectraWeights.h"
#include "AliAnalysisTaskMKBase.h"
#include "AliAnalysisTaskBaseWeights.h"

class AliAnalysisTaskBaseWeights;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskBaseWeights)
/// \endcond
//_____________________________________________________________________________

AliAnalysisTaskBaseWeights::AliAnalysisTaskBaseWeights()
    : AliAnalysisTaskMKBase()
    , fUseMCWeights(kTRUE)
    , fUseRandomSeed(kFALSE)
    , fRand(0)
    , fMCSpectraWeights(0)
    , fMCweight(1)
{
    // default contructor
}

//_____________________________________________________________________________

AliAnalysisTaskBaseWeights::AliAnalysisTaskBaseWeights(const char* name)
    : AliAnalysisTaskMKBase(name)
    , fUseMCWeights(kTRUE)
    , fUseRandomSeed(kFALSE)
    , fRand(0)
    , fMCSpectraWeights(0)
    , fMCweight(1)
{
    // constructor
}

//_____________________________________________________________________________

AliAnalysisTaskBaseWeights::~AliAnalysisTaskBaseWeights()
{
    if (fRand) { delete fRand; fRand=0; }
}

//_____________________________________________________________________________

void AliAnalysisTaskBaseWeights::BaseAddOutput()
{
    // add also all the default output from the base class
    AliAnalysisTaskMKBase::BaseAddOutput();

    // TODO add some control histograms here

    // if there are weights and they are used add them to the ouput
    if (fMCSpectraWeights && fUseMCWeights) {
        //TODO
        // removed after comment from patrick, only store the histogram (see below)
        //fOutputList->Add(fMCSpectraWeights);
        fOutputList->Add((TObject*)fMCSpectraWeights->GetHistMCGenPrimTrackParticles());
    }

}


//_____________________________________________________________________________

UInt_t AliAnalysisTaskBaseWeights::GetSeed()
{
    if (fUseRandomSeed) { return 0; }

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

void AliAnalysisTaskBaseWeights::LoopOverAllTracks(Int_t flag)
{
    // on data or if fUseMCWeights is not set we do not modify anything
    if (!fIsMC) { AliAnalysisTaskMKBase::LoopOverAllTracks(flag); }
    if (!fUseMCWeights) { AliAnalysisTaskMKBase::LoopOverAllTracks(flag); }

    fNTracksESD = fESD->GetNumberOfTracks();
    for (Int_t i = 0; i < fNTracksESD; i++) {
        fESDTrack = dynamic_cast<AliESDtrack*>(fESD->GetTrack(i));
        if (!fESDTrack) { Err("noESDtrack"); continue; }
        InitTrack();
        // get the scaling factor
        fMCweight = MCScalingFactor();
        Double_t s = fMCweight;
        while (s >= 1) {
            BaseAnaTrack(flag);
            s--;
        }
        if (s > 0) {
            if (!fRand) { fRand = new TRandom3(); }
            fRand->SetSeed(GetSeed());
            if (fRand->Rndm() < s) { BaseAnaTrack(flag); }
        }

    }
}


//_____________________________________________________________________________

void AliAnalysisTaskBaseWeights::LoopOverAllParticles(Int_t flag)
{
    // this method should not be called on data
    if (!fIsMC) return;
    // only if fUseMCWeights is set we do re weighting
    if (!fUseMCWeights) { AliAnalysisTaskMKBase::LoopOverAllParticles(flag); }

    fMCnTracks = fMC->GetNumberOfTracks();
    for (Int_t i = 0; i < fMCnTracks; i++) {
        fMCParticle  = dynamic_cast<AliMCParticle*>(fMC->GetTrack(i));
        if (!fMCParticle) { Err("noMCParticle"); continue; }
        fMCLabel = i;
        InitMCParticle();
        // get the scaling factor
        fMCweight = MCScalingFactor();
        Double_t s = fMCweight;
        while (s >= 1) {
            BaseAnaParticleMC(flag);
            s--;
        }
        if (s > 0) {
            if (!fRand) { fRand = new TRandom3(); }
            fRand->SetSeed(GetSeed());
            if (fRand->Rndm() < s) { BaseAnaParticleMC(flag); }
        }
    }
}

//_____________________________________________________________________________

Double_t AliAnalysisTaskBaseWeights::MCScalingFactor()
{
    // determine the MC scaling factor for the current particle

    // in case the weights are not activated, they are always unity
    if (!fUseMCWeights) { return 1.0; }

    // in case mcspectraweights are there we use them for primary particles
    if (fMCSpectraWeights && fMCisPrim) {
        // for now I pass the V0M multiplicity percentile
        // TODO check if this the correct one or what should be passed

        // TODO maybe ideal is to pass simply the event and let the SpectraWeights decide on which mult to use
        cout<<"fMCSpectraWeights->GetMCSpectraWeight(fMCParticle->Particle(), fMC)= "<<fMCSpectraWeights->GetMCSpectraWeight(fMCParticle->Particle(), fMC)<<endl; //DEBUG
        return fMCSpectraWeights->GetMCSpectraWeight(fMCParticle->Particle(), fMC);
    }

    // if there are no spectraweights or for secondaries we use the static tools function
    // TODO this has to be modified, currently values for the LCH17pq are hardcoded
    // TODO also the systematic varation should be set in the addtask somewhere
    // now always nominal values are used
    return AlidNdPtTools::MCScalingFactor(fMCParticle, fMC, 0);
}

//_____________________________________________________________________________

void AliAnalysisTaskBaseWeights::FillDefaultHistograms(Int_t step)
{
    // fill the spectra weights if needed
    // but do it before any event selection
    // because this might depent on the derived task

    //protection
    if (fMCSpectraWeights && step==0) {
        if(fMCSpectraWeights->GetTaskStatus() < AliMCSpectraWeights::TaskState::kMCSpectraObtained) {
            // for now I pass the V0M multiplicity percentile
            // TODO check if this the correct one or what should be passed
            // TODO maybe ideal is to pass simply the event and let the SpectraWeights decide on which mult to use
            fMCSpectraWeights->FillMCSpectra(fMC);
            cout<<"fMCSpectraWeights->FillMCSpectra(fMC);"<<endl; //DEBUG
        }
    }

    //also fill the default histograms
    AliAnalysisTaskMKBase::FillDefaultHistograms(step);

}

//_____________________________________________________________________________

AliAnalysisTaskBaseWeights* AliAnalysisTaskBaseWeights::AddTaskBaseWeights(const char* name, const char* outfile, const char* collisionSystem, Int_t sysFlag, const char* prevTrainOutputPath)
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskBaseWeights", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskBaseWeights", "This task requires an input event handler");
        return NULL;
    }

    // Setup output file
    //===========================================================================
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":";
    fileName += name;  // create a subfolder in the file
    if (outfile) { // if a finename is given, use that one
        fileName = TString(outfile);
    }


    // create the task
    //===========================================================================
    AliAnalysisTaskBaseWeights *task = new AliAnalysisTaskBaseWeights(name);
    if (!task) { return 0; }

    // configure the task
    //===========================================================================
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
    task->SetESDtrackCuts(0,AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));

    // configure the use of AliMCSpectraWeights
    //===========================================================================
    // collisionSystem is String "pp", "pPb", "XeXe", "PbPb"
    if (collisionSystem) {
        AliMCSpectraWeights* weights = new AliMCSpectraWeights(collisionSystem, "fMCSpectraWeights", (AliMCSpectraWeights::SysFlag)sysFlag);
        // root file with fHistMCGenPrimTrackParticle (THnF)
        // or used SetSavedListName
        // prevTrainOutputPath is string, for lego train activate the check box addtask to get alien connection for addtask
        // file is only needed for addtask
        //
        // file with fractions from data: expert input from Patrick:
        // /alice/cern.ch/user/p/phuhn/AllPublishedFractions.root
        // (path fixed in AliMCSpectraWeights code)
        //  used void SetDataFractionsFile(const char* file) to set local path

        if (prevTrainOutputPath) { weights->SetMCSpectraFile(prevTrainOutputPath); } // path to previous train output
        weights->Init();
        task->fMCSpectraWeights = weights;
    }


    // attach the task to the manager and configure in and ouput
    //===========================================================================
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

  return task;
}
