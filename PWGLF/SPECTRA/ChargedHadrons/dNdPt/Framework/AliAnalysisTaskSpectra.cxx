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
#include "AliAnalysisTaskSpectra.h"

class AliAnalysisTaskSpectra;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSpectra)
/// \endcond
//_____________________________________________________________________________

AliAnalysisTaskSpectra::AliAnalysisTaskSpectra() 
    : AliAnalysisTaskMKBase()
    , fHistEffCont(0)
    , fHistTrack(0)   
    , fHistEvent(0)
{
    // default contructor
 
}

//_____________________________________________________________________________

AliAnalysisTaskSpectra::AliAnalysisTaskSpectra(const char* name) 
    : AliAnalysisTaskMKBase(name)
    , fHistEffCont(0)
    , fHistTrack(0)
    , fHistEvent(0)
{
    // constructor

}

//_____________________________________________________________________________

AliAnalysisTaskSpectra::~AliAnalysisTaskSpectra()
{
    // destructor
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectra::AddOutput()
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
    AddAxis("nAcc","mult6kcoarse"); 
    AddAxis("pt");         
    AddAxis("Q",3,-1.5,1.5);        
    fHistTrack = CreateHist("fHistTrack");
    fOutputList->Add(fHistTrack);
        
    AddAxis("cent");    
    AddAxis("nAcc","mult6kfine");
    AddAxis("zV",8,-20,20);
    fHistEvent = CreateHist("fHistEvent");
    fOutputList->Add(fHistEvent);    
    
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskSpectra::IsEventSelected()
{
    return fIsAcceptedAliEventCuts;
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectra::AnaEvent()
{
   
   LoopOverAllTracks();
   if (fIsMC) LoopOverAllParticles();
   
   FillHist(fHistEvent, fMultPercentileV0M, fNTracksAcc, fZv);
   
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectra::AnaTrack(Int_t flag)
{
    if (!fAcceptTrackM) return;
    
    FillHist(fHistTrack, fMultPercentileV0M, fNTracksAcc, fPt, fChargeSign);
}

//_____________________________________________________________________________

void AliAnalysisTaskSpectra::AnaTrackMC(Int_t flag)
{
    if (!fAcceptTrackM) return;
     
    if (fMCParticleType==AlidNdPtTools::kOther) { Log("RecTrack.PDG.",fMCPDGCode); }
    if (TMath::Abs(fMCQ > 1)) { Log("RecTrack.Q>1.PDG.",fMCPDGCode); }
    
    FillHist(fHistEffCont, fMultPercentileV0M, fNTracksAcc, fMCPt, fMCChargeSign, fMCParticleType, fMCProdcutionType); 

}

//_____________________________________________________________________________

void AliAnalysisTaskSpectra::AnaParticleMC(Int_t flag)
{            
    if (!fMCisPrim) return;    
    if (!fMCIsCharged) return;    
    if (TMath::Abs(fMCEta) > 0.8) return;    
    
    if (fMCParticleType==AlidNdPtTools::kOther) { Log("GenPrim.PDG.",fMCPDGCode); }
    if (TMath::Abs(fMCQ > 1)) { Log("GenPrim.Q>1.PDG.",fMCPDGCode); }
    
    FillHist(fHistEffCont, fMultPercentileV0M, fNTracksAcc, fMCPt, fMCChargeSign, fMCParticleType, 3);     

}

//_____________________________________________________________________________

AliAnalysisTaskSpectra* AliAnalysisTaskSpectra::AddTaskSpectra(const char* name, const char* outfile) 
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskSpectra", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskSpectra", "This task requires an input event handler");
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
    AliAnalysisTaskSpectra *task = new AliAnalysisTaskSpectra(name);  
    if (!task) { return 0; }
    
    // configure the task
    //===========================================================================
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);    
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
//     task->SetESDtrackCuts(0,AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
    task->SetNeedEventMult(kTRUE);
    task->SetNeedEventCent(kTRUE);

    // attach the task to the manager and configure in and ouput
    //===========================================================================
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
  return task;
}
