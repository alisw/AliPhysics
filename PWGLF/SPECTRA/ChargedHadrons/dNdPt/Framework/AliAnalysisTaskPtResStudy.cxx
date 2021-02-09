#include <iostream>
#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
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
#include "AliAnalysisTaskPtResStudy.h"

class AliAnalysisTaskPtResStudy;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskPtResStudy)
/// \endcond
//_____________________________________________________________________________

AliAnalysisTaskPtResStudy::AliAnalysisTaskPtResStudy() 
    : AliAnalysisTaskMKBase()
    , fHistPtResCov(0)
    , fHistPtResCovHighPt(0)
    , fHistPtResMC(0)
    , fHistPtRes(0)
{
    // default contructor

}

//_____________________________________________________________________________

AliAnalysisTaskPtResStudy::AliAnalysisTaskPtResStudy(const char* name) 
    : AliAnalysisTaskMKBase(name)
    , fHistPtResCov(0)  
    , fHistPtResCovHighPt(0)
    , fHistPtResMC(0)
    , fHistPtRes(0)
{
    // constructor

}

//_____________________________________________________________________________

AliAnalysisTaskPtResStudy::~AliAnalysisTaskPtResStudy()
{
    // destructor
}

//_____________________________________________________________________________

void AliAnalysisTaskPtResStudy::AddOutput()
{    
    // data histogram for covariance matrix entries at high pt
    // signed 1/pt : sigma(1/pt) : ntracks
    AddAxis("signed1pt","1/pT",2000,-10,10);    
    AddAxis("sigma1pt","#sigma(1/pT)",1000,0,0.1);    
    AddAxis("nTracks","mult6kcoarse");
    fHistPtResCov = CreateHist("fHistPtResCov");
    fOutputList->Add(fHistPtResCov);

    // data histogram for covariance matrix entries at high pt
    // 1/pt : sigma(1/pt) : ntracks
    AddAxis("1pt","1/pT",200,0,0.2);    
    AddAxis("sigma1pt","#sigma(1/pT)",2000,0,0.02);    
    AddAxis("nTracks","mult6kcoarse");
    fHistPtResCovHighPt = CreateHist("fHistPtResCovHighPt");
    fOutputList->Add(fHistPtResCovHighPt);
    
    
    // pt response in small pT bins for low pt
    // pt : ptMC : nTracks
    AddAxis("pt",1000,0,10);    
    AddAxis("ptMC",1000,0,10);    
    AddAxis("nTracks","mult6kcoarse");
    fHistPtResMC = CreateHist("fHistPtResMC");
    fOutputList->Add(fHistPtResMC);
    
    // mc histogram for covariance matrix entries
    // pt : ptmc : ntracks
    AddAxis("pt","pt");    
    AddAxis("ptMC","ptMC");
    AddAxis("sigmapt","#sigma(1/pT)*pT",100,0,0.1);
    AddAxis("deltapt","(ptmc/pt-1)",200,-0.1,0.1);
    AddAxis("nTracks","mult6kcoarse");
    fHistPtRes = CreateHist("fHistPtRes");
    fOutputList->Add(fHistPtRes);
    
    
}


//_____________________________________________________________________________

Bool_t AliAnalysisTaskPtResStudy::IsEventSelected()
{
    return fIsAcceptedAliEventCuts;
}

//_____________________________________________________________________________

void AliAnalysisTaskPtResStudy::AnaEvent()
{
    LoopOverAllTracks();
}

//_____________________________________________________________________________

void AliAnalysisTaskPtResStudy::AnaTrack(Int_t flag)
{
    if (!fAcceptTrack[0]) return;
    if (f1Pt<1.0)  { FillHist(fHistPtResCov,       fSigned1Pt, fSigma1Pt, fNTracksAcc); }
    if (f1Pt<0.2)  { FillHist(fHistPtResCovHighPt, f1Pt, fSigma1Pt, fNTracksAcc); }
    if (fMCPt<10.) { FillHist(fHistPtResMC,        fPt, fMCPt, fNTracksAcc); }
    FillHist(fHistPtRes,          fPt, fMCPt, fSigma1Pt*fPt, fMCPt/fPt-1, fNTracksAcc);
}

//_____________________________________________________________________________

AliAnalysisTaskPtResStudy* AliAnalysisTaskPtResStudy::AddTaskPtResStudy(const char* name, const char* outfile) 
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskPtResStudy", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskPtResStudy", "This task requires an input event handler");
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
    AliAnalysisTaskPtResStudy *task = new AliAnalysisTaskPtResStudy(name);  
    if (!task) { return 0; }
    
    // configure the task
    //===========================================================================
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);    
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));    
    task->SetESDtrackCuts(0,AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));        
    task->SetNeedEventMult(kTRUE);
    task->SetNeedTrackTPC(kTRUE);

    // attach the task to the manager and configure in and ouput
    //===========================================================================
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
  return task;
}
