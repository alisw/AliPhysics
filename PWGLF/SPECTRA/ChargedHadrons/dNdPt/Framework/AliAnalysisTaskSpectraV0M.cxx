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
#include "AliAnalysisTaskSpectraV0M.h"

class AliAnalysisTaskSpectraV0M;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSpectraV0M)
/// \endcond
//_____________________________________________________________________________

AliAnalysisTaskSpectraV0M::AliAnalysisTaskSpectraV0M() 
    : AliAnalysisTaskMKBase()
    , fHistEffCont(0)
    , fHistTrack(0)   
    , fHistEvent(0)
{
    // default contructor
}

//_____________________________________________________________________________

AliAnalysisTaskSpectraV0M::AliAnalysisTaskSpectraV0M(const char* name) 
    : AliAnalysisTaskMKBase(name)
    , fHistEffCont(0)
    , fHistTrack(0)
    , fHistEvent(0)
{
    // constructor
}

//_____________________________________________________________________________

AliAnalysisTaskSpectraV0M::~AliAnalysisTaskSpectraV0M()
{
    // destructor
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraV0M::AddOutput()
{        
    AddAxis("cent");
    AddAxis("nAcc","mult6kcoarse");
    AddAxis("MCpT","pt");
    AddAxis("MCQ",3,-1.5,1.5);        
    AddAxis("MCpid",10,-0.5,9.5);  // 0=e, 1=mu, 2=pi, 3=K, 4=p, 6=sigmaP, 7=sigmaM, 8=xi, 9=omega, 5=other
    AddAxis("MCinfo",4,-0.5,3.5);  // 0=prim, 1=decay 2=material, 3=genprim    
    fHistEffCont = CreateHist("fHistEffCont");
    fOutputList->Add(fHistEffCont);
    
    AddAxis("cent");
    AddAxis("multV0","mult6kfine");
    AddAxis("pt");             
    fHistTrack = CreateHist("fHistTrack");
    fOutputList->Add(fHistTrack);
        
    
    AddAxis("cent");
    AddAxis("multV0","mult6kfine");
    AddAxis("nAcc","mult6kcoarse");
    fHistEvent = CreateHist("fHistEvent");
    fOutputList->Add(fHistEvent);    
    
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskSpectraV0M::IsEventSelected()
{
    return fIsAcceptedAliEventCuts;
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraV0M::AnaEvent()
{
   
   LoopOverAllTracks();
//    if (fIsMC) LoopOverAllParticles();
   
   FillHist(fHistEvent, fMultPercentileV0M, fMultV0MmultSelection, fNTracksAcc);
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraV0M::AnaTrack(Int_t flag)
{
    if (!fAcceptTrackM) return;
    
    FillHist(fHistTrack, fMultPercentileV0M, fMultV0MmultSelection, fPt);
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraV0M::AnaTrackMC(Int_t flag)
{
    if (!fAcceptTrackM) return;
     
    if (fMCParticleType==AlidNdPtTools::kOther) { Log("RecTrack.PDG.",fMCPDGCode); }
    if (TMath::Abs(fMCQ > 1)) { Log("RecTrack.Q>1.PDG.",fMCPDGCode); }
    
    FillHist(fHistEffCont, fMultPercentileV0M, fNTracksAcc, fMCPt, fMCChargeSign, fMCParticleType, fMCProdcutionType); 

}

//_____________________________________________________________________________

void AliAnalysisTaskSpectraV0M::AnaParticleMC(Int_t flag)
{            
    if (!fMCisPrim) return;    
    if (!fMCIsCharged) return;    
    if (TMath::Abs(fMCEta) > 0.8) return;    
    
    if (fMCParticleType==AlidNdPtTools::kOther) { Log("GenPrim.PDG.",fMCPDGCode); }
    if (TMath::Abs(fMCQ > 1)) { Log("GenPrim.Q>1.PDG.",fMCPDGCode); }
    
    FillHist(fHistEffCont, fMultPercentileV0M, fNTracksAcc, fMCPt, fMCChargeSign, fMCParticleType, 3);     

}

//_____________________________________________________________________________

AliAnalysisTaskSpectraV0M* AliAnalysisTaskSpectraV0M::AddTaskSpectraV0M(const char* name, const char* outfile) 
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskSpectraV0M", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskSpectraV0M", "This task requires an input event handler");
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
    AliAnalysisTaskSpectraV0M *task = new AliAnalysisTaskSpectraV0M(name);  
    if (!task) { return 0; }
    
    // configure the task
    //===========================================================================
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);    
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
//     task->SetESDtrackCuts(0,AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
    
    // attach the task to the manager and configure in and ouput
    //===========================================================================
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
  return task;
}
