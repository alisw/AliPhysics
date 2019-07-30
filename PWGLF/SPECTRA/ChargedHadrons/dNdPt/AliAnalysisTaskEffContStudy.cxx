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
#include "AliAnalysisTaskEffContStudy.h"

class AliAnalysisTaskEffContStudy;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEffContStudy)
/// \endcond
//_____________________________________________________________________________

AliAnalysisTaskEffContStudy::AliAnalysisTaskEffContStudy() 
    : AliAnalysisTaskMKBase()
    , fHistEffCont(0)
    , fHistEffContScaled(0)   
{
    // default contructor
}

//_____________________________________________________________________________

AliAnalysisTaskEffContStudy::AliAnalysisTaskEffContStudy(const char* name) 
    : AliAnalysisTaskMKBase(name)
    , fHistEffCont(0)
    , fHistEffContScaled(0)
{
    // constructor
}

//_____________________________________________________________________________

AliAnalysisTaskEffContStudy::~AliAnalysisTaskEffContStudy()
{
    // destructor
}

//_____________________________________________________________________________

void AliAnalysisTaskEffContStudy::AddOutput()
{        
    AddAxis("MCpT","pt");
    AddAxis("MCpid",10,-0.5,9.5);  // 0=e, 1=mu, 2=pi, 3=K, 4=p, 6=sigmaP, 7=sigmaM, 8=xi, 9=omega, 5=other
    AddAxis("MCinfo",4,-0.5,3.5); // 0=prim, 1=decay 2=material, 3=genprim
    AddAxis("MCQ",3,-1.5,1.5);
    AddAxis("nTracks","mult6kcoarse");    
    AddAxis("nPrimMCpT150","mult6kcoarse");   
    AddAxis("nPrimMCallpT","mult6kcoarse");   
    fHistEffCont = CreateHist("fHistEffCont");
    fOutputList->Add(fHistEffCont);
    
    AddAxis("MCpT","pt");
    AddAxis("MCpid",10,-0.5,9.5);  // 0=e, 1=mu, 2=pi, 3=K, 4=p, 6=sigmaP, 7=sigmaM, 8=xi, 9=omega, 5=other
    AddAxis("MCinfo",4,-0.5,3.5); // 0=prim, 1=decay 2=material, 3=genprim
    AddAxis("MCQ",3,-1.5,1.5);
    AddAxis("nTracks","mult6kcoarse");
    fHistEffContScaled = CreateHist("fHistEffContScaled");
    fOutputList->Add(fHistEffContScaled);
    
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskEffContStudy::IsEventSelected()
{
    return fIsAcceptedAliEventCuts;
}

//_____________________________________________________________________________

void AliAnalysisTaskEffContStudy::AnaEvent()
{
   
   LoopOverAllTracks();
   if (fIsMC) LoopOverAllParticles();
   
}

//_____________________________________________________________________________

void AliAnalysisTaskEffContStudy::AnaTrackMC(Int_t flag)
{
    if (!fAcceptTrack[0]) return;
    FillHist(fHistEffCont, fMCPt, fMCParticleType, fMCProdcutionType, fMCChargeSign, fNTracksAcc, fMCnPrimPtCut, fMCnPrim08); 
    
    if (fMCParticleType==AlidNdPtTools::kOther) { Log("RecTrack.PDG.",fMCPDGCode); }
    if (TMath::Abs(fMCQ > 1)) { Log("RecTrack.Q>1.PDG.",fMCPDGCode); }    
    
    Double_t s = AlidNdPtTools::MCScalingFactor(fMCProdcutionType,fMCParticleType, fMCPt); 
    
    while (s >= 1) {
        FillHist(fHistEffContScaled, fMCPt, fMCParticleType, fMCProdcutionType, fMCChargeSign, fNTracksAcc); 
        s--;
    }
    if (s > 0) {
        if (gRandom->Rndm() < s) { FillHist(fHistEffContScaled, fMCPt, fMCParticleType, fMCProdcutionType, fMCChargeSign, fNTracksAcc); }
    }
}

//_____________________________________________________________________________

void AliAnalysisTaskEffContStudy::AnaParticleMC(Int_t flag)
{            
    if (!fMCisPrim) return;    
    if (!fMCIsCharged) return;    
    if (TMath::Abs(fMCEta) > 0.8) return;    
    
    if (fMCParticleType==AlidNdPtTools::kOther) { Log("GenPrim.PDG.",fMCPDGCode); }
    if (TMath::Abs(fMCQ > 1)) { Log("GenPrim.Q>1.PDG.",fMCPDGCode); }
    
    FillHist(fHistEffCont, fMCPt, fMCParticleType, 3, fMCChargeSign, fNTracksAcc, fMCnPrimPtCut, fMCnPrim08);        
    
    Double_t s = AlidNdPtTools::MCScalingFactor(fMCProdcutionType,fMCParticleType, fMCPt);    
        
    while (s >= 1) {
        FillHist(fHistEffContScaled, fMCPt, fMCParticleType, 3, fMCChargeSign, fNTracksAcc); 
        s--;
    }
    if (s > 0) {
        if (gRandom->Rndm() < s) { FillHist(fHistEffContScaled, fMCPt, fMCParticleType, 3, fMCChargeSign, fNTracksAcc); }
    }
}

//_____________________________________________________________________________

AliAnalysisTaskEffContStudy* AliAnalysisTaskEffContStudy::AddTaskEffContStudy(const char* name, const char* outfile) 
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskEffContStudy", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskEffContStudy", "This task requires an input event handler");
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
    AliAnalysisTaskEffContStudy *task = new AliAnalysisTaskEffContStudy(name);  
    if (!task) { return 0; }
    
    // configure the task
    //===========================================================================
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);    
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
    task->SetESDtrackCuts(0,AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
    
    // attach the task to the manager and configure in and ouput
    //===========================================================================
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
  return task;
}