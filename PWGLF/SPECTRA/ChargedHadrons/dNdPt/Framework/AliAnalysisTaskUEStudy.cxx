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
#include "AliAnalysisTaskUEStudy.h"

class AliAnalysisTaskUEStudy;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskUEStudy)
/// \endcond
//_____________________________________________________________________________

AliAnalysisTaskUEStudy::AliAnalysisTaskUEStudy() 
    : AliAnalysisTaskMKBase()
    , fPtMax(-1)
    , fPtMaxPhi(-1)
    , fMCPtMax(-1)
    , fMCPtMaxPhi(-1)
    , fNTracksTowards(-1)
    , fNTracksTransverse(-1)
    , fNTracksAway(-1)
    , fMCNChTowards(-1)
    , fMCNChTransverse(-1)
    , fMCNChAway(-1)        
    , fHistUETracks(0)
    , fHistUE(0)   
    , fHistUETracksMC(0)
    , fHistUEMC(0)   
    , fHistUEPhiRes(0)
{
    // default contructor

}

//_____________________________________________________________________________

AliAnalysisTaskUEStudy::AliAnalysisTaskUEStudy(const char* name) 
    : AliAnalysisTaskMKBase(name)
    , fPtMax(-1)
    , fPtMaxPhi(-1)
    , fMCPtMax(-1)
    , fMCPtMaxPhi(-1)
    , fNTracksTowards(-1)
    , fNTracksTransverse(-1)
    , fNTracksAway(-1)
    , fMCNChTowards(-1)
    , fMCNChTransverse(-1)
    , fMCNChAway(-1)        
    , fHistUETracks(0)
    , fHistUE(0)   
    , fHistUETracksMC(0)
    , fHistUEMC(0)   
    , fHistUEPhiRes(0)
{
    // constructor
 
}

//_____________________________________________________________________________

AliAnalysisTaskUEStudy::~AliAnalysisTaskUEStudy()
{
    // destructor
}

//_____________________________________________________________________________

void AliAnalysisTaskUEStudy::AddOutput()
{   
    //pt:ptmax:ntracksaway:cosdphi:ntracksaway:ntrackstransverse:ntrackstowards
    AddAxis("pT","pt");
    AddAxis("pTmax","ptveryfew");
    AddAxis("cos(deltaphi)",3,-1.5,1.5); // 3bins: 1=away, 2=transverse, 3=towards    
    AddAxis("nTracksAway","mult6kfine");    
    AddAxis("nTracksTransverse","mult6kfine");
    AddAxis("nTracksTowards","mult6kfine");
    fHistUETracks = CreateHist("fHistUETracks");
    fOutputList->Add(fHistUETracks);
    
    //pt:ptmax:ntracksaway:cosdphi:ntracksaway:ntrackstransverse:ntrackstowards
    AddAxis("pTMC","pt");    
    AddAxis("pTmaxMC","ptveryfew");
    AddAxis("cos(deltaphiMC)",3,-1.5,1.5); // 3bins: 1=away, 2=transverse, 3=towards    
    AddAxis("multAway","mult6kfine");    
    AddAxis("multTansverse","mult6kfine");
    AddAxis("multTowards","mult6kfine");
    fHistUETracksMC = CreateHist("fHistUETracksMC");
    fOutputList->Add(fHistUETracksMC);
    
    //ptmax:ntracksaway;ntrackstransverse:ntracktowards:ntracks:mcptmax:nchawa:nchtransverse:nchtowards:nch
    AddAxis("pTmax","ptveryfew");    
    AddAxis("nTracksAway","mult6kfine");    
    AddAxis("nTracksTransverse","mult6kfine");
    AddAxis("nTracksTowards","mult6kfine");
//     AddAxis("nTracks","mult6kfine");
    AddAxis("pTmaxMC","ptveryfew");    
    AddAxis("nTracksAwayMC","mult6kfine");    
    AddAxis("nTracksTransverseMC","mult6kfine");
    AddAxis("nTracksTowardsMC","mult6kfine");
//     AddAxis("nTracksMC","mult6kfine");
    fHistUE = CreateHist("fHistUE");
    fOutputList->Add(fHistUE);    
    
    //ptmax:ptmaxmc:phimax:phimaxmc:deltaphi    
    AddAxis("pTmax","ptveryfew");    
    AddAxis("pTmaxMC","ptveryfew");            
    AddAxis("PhiMax",180,0,2*TMath::Pi());
    AddAxis("PhiMaxMC",180,0,2*TMath::Pi());
    AddAxis("DeltaPhi",360,-0.1*TMath::Pi(),0.1*TMath::Pi());
    fHistUEPhiRes = CreateHist("fHistUEPhiRes");
    fOutputList->Add(fHistUEPhiRes);
}


//_____________________________________________________________________________

Bool_t AliAnalysisTaskUEStudy::IsEventSelected()
{
    return fIsAcceptedAliEventCuts;
}

//_____________________________________________________________________________

void AliAnalysisTaskUEStudy::AnaEvent()
{      
   
    fPtMax = -1;
    fMCPtMax = -1;    
    fNTracksTowards = 0;
    fNTracksTransverse = 0;
    fNTracksAway = 0;
    fMCNChTowards = 0;
    fMCNChTransverse = 0;
    fMCNChAway = 0;
   
    // 1st loop to find highest pt particle    
    LoopOverAllTracks(1);
    LoopOverAllParticles(1);    
    // 2nd loop to cout multiplicities
    LoopOverAllTracks(2);
    LoopOverAllParticles(2);    
    // 3nd loop to fill histograms
    LoopOverAllTracks(3);
    LoopOverAllParticles(3);
    
    // fill event based histograms
    
    //ptmax:ntracksaway;ntrackstransverse:ntracktowards:ntracks:mcptmax:nchawa:nchtransverse:nchtowards:nch
    FillHist(fHistUE,fPtMax,fNTracksAway,fNTracksTransverse,fNTracksTowards,fMCPtMax,fMCNChAway,fMCNChTransverse,fMCNChTowards);
   
    //ptmax:ptmaxmc:phimax:phimaxmc:deltaphi    
    Double_t deltaphi = fPtMaxPhi-fMCPtMaxPhi;
    AlidNdPtTools::Range1Pi(deltaphi);
   
    FillHist(fHistUEPhiRes,fPtMax,fMCPtMax,fPtMaxPhi,fMCPtMaxPhi,deltaphi);
   
   
}

//_____________________________________________________________________________

void AliAnalysisTaskUEStudy::AnaTrack(Int_t flag)
{
    if (!fAcceptTrack[0]) return;    
        
    if (flag == 1) {             
        if (fPt>fPtMax) { fPtMax = fPt; fPtMaxPhi = fPhi; }
        return;
    }
    
    Double_t cosdphi = TMath::Cos(fPhi-fPtMaxPhi);
    if (flag == 2) { 
        if (cosdphi < -0.5) {
            fNTracksAway++;
        } else if (cosdphi >= 0.5) {
            fNTracksTowards++;
        } else {
            fNTracksTransverse++;
        }
    }
    
    if (flag == 3) {
        //pt:ptmax:ntracksaway:cosdphi:ntracksaway:ntrackstransverse:ntrackstowards
        FillHist(fHistUETracks,fPt,fPtMax,cosdphi,fNTracksAway,fNTracksTransverse,fNTracksTowards);        
    }
}

//_____________________________________________________________________________

void AliAnalysisTaskUEStudy::AnaParticleMC(Int_t flag)
{
    if (!fMCisPrim) return;    
    if (!fMCIsCharged) return;    
    if (TMath::Abs(fMCEta) > 0.8) return;  
    
        
    if (flag == 1) {        
        if (fMCPt>fMCPtMax) { fMCPtMax = fMCPt; fMCPtMaxPhi = fMCPhi; }
        return;
    }
    
    Double_t cosdphi = TMath::Cos(fMCPhi-fMCPtMaxPhi);
    if (flag == 2) { 
        if (cosdphi < -0.5) {
            fMCNChAway++;
        } else if (cosdphi >= 0.5) {
            fMCNChTransverse++;
        } else {
            fMCNChTransverse++;
        }
    }
    
    if (flag == 3) {
        //pt:ptmax:ntracksaway:cosdphi:ntracksaway:ntrackstransverse:ntrackstowards
        FillHist(fHistUETracksMC,fMCPt,fMCPtMax,cosdphi,fMCNChAway,fMCNChTransverse,fMCNChTowards);        
    }
}

//_____________________________________________________________________________

AliAnalysisTaskUEStudy* AliAnalysisTaskUEStudy::AddTaskUEStudy(const char* name, const char* outfile) 
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskUEStudy", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskUEStudy", "This task requires an input event handler");
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
    AliAnalysisTaskUEStudy *task = new AliAnalysisTaskUEStudy(name);  
    if (!task) { return 0; }
    
    // configure the task
    //===========================================================================
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);    
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
    task->SetESDtrackCuts(0,AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
    task->SetNeedEventMult(kTRUE);
    task->SetNeedEventCent(kTRUE);
    task->SetNeedTrackTPC(kTRUE);
    // attach the task to the manager and configure in and ouput
    //===========================================================================
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
  return task;
}
