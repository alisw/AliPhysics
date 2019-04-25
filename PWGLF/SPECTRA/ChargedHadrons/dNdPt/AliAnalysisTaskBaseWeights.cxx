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
        // the implementation in AlidNdPtTools is temporary for testing puposes
        fMCweight = AlidNdPtTools::MCScalingFactor(fMCProdcutionType,fMCParticleType, fMCPt);
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
        // the implementation in AlidNdPtTools is temporary for testing puposes
        fMCweight = AlidNdPtTools::MCScalingFactor(fMCProdcutionType,fMCParticleType, fMCPt);
        Double_t s = fMCweight;
        while (s >= 1) {
            BaseAnaParticleMC(flag);
            s--;
        }
        if (s > 0) {
            fRand->SetSeed(GetSeed());
            if (fRand->Rndm() < s) { BaseAnaParticleMC(flag); }
        }        
    }
}

//_____________________________________________________________________________

AliAnalysisTaskBaseWeights* AliAnalysisTaskBaseWeights::AddTaskBaseWeights(const char* name, const char* outfile) 
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
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("default"));
    task->SetESDtrackCuts(0,AlidNdPtTools::CreateESDtrackCuts("default"));
    
    // attach the task to the manager and configure in and ouput
    //===========================================================================
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
  return task;
}