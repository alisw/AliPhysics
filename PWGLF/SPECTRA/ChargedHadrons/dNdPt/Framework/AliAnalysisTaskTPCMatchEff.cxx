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
#include "AliAnalysisTaskTPCMatchEff.h"

class AliAnalysisTaskTPCMatchEff;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskTPCMatchEff)
/// \endcond
//_____________________________________________________________________________

AliAnalysisTaskTPCMatchEff::AliAnalysisTaskTPCMatchEff() 
    : AliAnalysisTaskMKBase()
    , fHistMCMatchEff(0)
    , fHistDATAMatchEff(0)
{
    // default contructor
}

//_____________________________________________________________________________

AliAnalysisTaskTPCMatchEff::AliAnalysisTaskTPCMatchEff(const char* name) 
    : AliAnalysisTaskMKBase(name)
    , fHistMCMatchEff(0)
    , fHistDATAMatchEff(0)   
{
    // constructor
}

//_____________________________________________________________________________

AliAnalysisTaskTPCMatchEff::~AliAnalysisTaskTPCMatchEff()
{
    // destructor
}

//_____________________________________________________________________________

void AliAnalysisTaskTPCMatchEff::AddOutput()
{    
    //mult:trackpt:tpcpt:tpceta:tpcphi:mcpt:sec:cut1:cut2:cut3:cut4
    AddAxis("centrality","V0M centrality","cent");
    AddAxis("trackpt","p_{t,track} (GeV/c)","pt");
    AddAxis("tpcpt","p_{t,TPC} (GeV/c)","pt");    
    AddAxis("tpceta","#eta_{TPC}",10,-1.,+1.);
    AddAxis("tpcphi","#phi_{TPC}",18,0.,2*TMath::Pi());
    if (fESDtrackCuts[1]) AddAxis("cut1",2,-0.5,1.5);
    if (fESDtrackCuts[2]) AddAxis("cut2",2,-0.5,1.5);
    if (fESDtrackCuts[3]) AddAxis("cut3",2,-0.5,1.5);  
    if (fESDtrackCuts[4]) AddAxis("cut4",2,-0.5,1.5);  
    AddAxis("MCinfo",3,-0.5,2.5); // 0=prim, 1=decay 2=material       
    AddAxis("mcpt","p_{t,MC} (GeV/c)","pt");
    fHistMCMatchEff = CreateHist("fHistMCMatchEff");
    fOutputList->Add(fHistMCMatchEff);  
    
    //cent:trackpt:tpcpt:tpceta:tpcphi:cut1:cut2:cut3:cut4
    AddAxis("centrality","V0M centrality","cent");
    AddAxis("trackpt","p_{t,track} (GeV/c)","pt");    
    AddAxis("tpcpt","p_{t,TPC} (GeV/c)","pt");    
    AddAxis("tpceta","#eta_{MC}",10,-1.,+1.);
    AddAxis("tpcphi","#phi_{MC}",18,0.,2*TMath::Pi());    
    if (fESDtrackCuts[1]) AddAxis("cut1",2,-0.5,1.5);
    if (fESDtrackCuts[2]) AddAxis("cut2",2,-0.5,1.5);
    if (fESDtrackCuts[3]) AddAxis("cut3",2,-0.5,1.5);  
    if (fESDtrackCuts[4]) AddAxis("cut4",2,-0.5,1.5);  
    fHistDATAMatchEff = CreateHist("fHistDATAMatchEff");
    fOutputList->Add(fHistDATAMatchEff);  
    
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskTPCMatchEff::IsEventSelected()
{
    return fIsAcceptedAliEventCuts;
}

//_____________________________________________________________________________

void AliAnalysisTaskTPCMatchEff::AnaEvent()
{
   LoopOverAllTracks();
}

//_____________________________________________________________________________

void AliAnalysisTaskTPCMatchEff::AnaTrackMC(Int_t flag)
{     
    if (!fAcceptTrack[0]) return;    
    FillHist(fHistMCMatchEff, fMultPercentileV0M, fPt, fPtInnerTPC, fEtaInnerTPC, fPhiInnerTPC, fAcceptTrack[1], fAcceptTrack[2], fAcceptTrack[3], fAcceptTrack[4], fMCPrimSec, fMCPt);
}

//_____________________________________________________________________________

void AliAnalysisTaskTPCMatchEff::AnaTrackDATA(Int_t flag)
{     
    if (!fAcceptTrack[0]) return;    
    FillHist(fHistDATAMatchEff, fMultPercentileV0M, fPt, fPtInnerTPC, fEtaInnerTPC, fPhiInnerTPC, fAcceptTrack[1], fAcceptTrack[2], fAcceptTrack[3], fAcceptTrack[4]);
}

//_____________________________________________________________________________

AliAnalysisTaskTPCMatchEff* AliAnalysisTaskTPCMatchEff::AddTaskTPCMatchEff(const char* name, const char* outfile) 
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskTPCMatchEff", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskTPCMatchEff", "This task requires an input event handler");
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
    AliAnalysisTaskTPCMatchEff *task = new AliAnalysisTaskTPCMatchEff(name);  
    if (!task) { return 0; }
    
    // configure the task
    //===========================================================================
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);    
    //task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
    task->SetESDtrackCuts(0,AlidNdPtTools::CreateESDtrackCuts("TPCgeoEta08"));    
    task->SetESDtrackCuts(1,AlidNdPtTools::CreateESDtrackCuts("TPCgeo+ITShitEta08"));
    task->SetESDtrackCuts(2,AlidNdPtTools::CreateESDtrackCuts("TPCgeo+ITSrefitEta08"));
    task->SetESDtrackCuts(3,AlidNdPtTools::CreateESDtrackCuts("TPCgeo+SPDhitEta08"));
    task->SetESDtrackCuts(4,AlidNdPtTools::CreateESDtrackCuts("TPCgeo+ITSrefit+SPDhitEta08"));
    
    
    // attach the task to the manager and configure in and ouput
    //===========================================================================
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
  return task;
}