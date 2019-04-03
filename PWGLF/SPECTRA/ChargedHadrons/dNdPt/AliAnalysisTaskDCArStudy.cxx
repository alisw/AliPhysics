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
#include "AliAnalysisTaskDCArStudy.h"

class AliAnalysisTaskDCArStudy;

using namespace std;

ClassImp(AliAnalysisTaskDCArStudy)

//_____________________________________________________________________________

AliAnalysisTaskDCArStudy::AliAnalysisTaskDCArStudy() 
    : AliAnalysisTaskMKBase()
    , fHistDCA(0)
    , fHistDCATPC(0)
{
    // default contructor
}

//_____________________________________________________________________________

AliAnalysisTaskDCArStudy::AliAnalysisTaskDCArStudy(const char* name) 
    : AliAnalysisTaskMKBase(name)
    , fHistDCA(0)    
    , fHistDCATPC(0)
{
    // constructor
}

//_____________________________________________________________________________

AliAnalysisTaskDCArStudy::~AliAnalysisTaskDCArStudy()
{
    // destructor
}

//_____________________________________________________________________________

void AliAnalysisTaskDCArStudy::AddOutput()
{    
    //dcar:pt:mult:mcinfo
    AddAxis("DCAxy",10000,-1,1);
    AddAxis("pt");    
    AddAxis("nTracks","mult6kfine");
    AddAxis("MCinfo",3,-0.5,2.5); // 0=prim, 1=decay 2=material
    fHistDCA = CreateHist("fHistDCA");
    fOutputList->Add(fHistDCA);
    
    //dcar:pt:mult:mcinfo
    AddAxis("DCAxy",10000,-20,20); 
    AddAxis("TPCpt","pt");    
    AddAxis("nTracks","mult6kfine");
    AddAxis("MCinfo",3,-0.5,2.5); // 0=prim, 1=decay 2=material
    fHistDCATPC = CreateHist("fHistDCATPC");
    fOutputList->Add(fHistDCATPC);
    
     //dcar:pt:mult:mcinfo
  
}

//_____________________________________________________________________________

void AliAnalysisTaskDCArStudy::AnaEvent()
{
   InitEvent();
   InitEventMult();
   InitEventCent();
   InitMCEvent();
   LoopOverAllTracks();
   
}

//_____________________________________________________________________________

void AliAnalysisTaskDCArStudy::AnaTrack()
{    
    InitTrack();
    InitMCTrack();
    InitTrackIP();
    InitTrackTPC();
    if (fESDtrackCuts[0]->AcceptTrack(fESDTrack)) { FillHist(fHistDCATPC, fDCArTPC, fPtInnerTPC, fNTracksAcc, fMCPrimSec); }
    if (fESDtrackCuts[1]->AcceptTrack(fESDTrack)) { FillHist(fHistDCA, fDCAr, fPt, fNTracksAcc, fMCPrimSec); }
}

//_____________________________________________________________________________

AliAnalysisTaskDCArStudy* AliAnalysisTaskDCArStudy::AddTaskDCArStudy(const char* name) 
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskDCArStudy", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskDCArStudy", "This task requires an input event handler");
        return NULL;
    }
    
    // Setup output file
    //===========================================================================
    TString fileName = AliAnalysisManager::GetCommonFileName();
    //fileName += ":TaskDCArStudy";      // create a subfolder in the file
    fileName = TString("out_");
    fileName += name;
    fileName += ".root";
    //if (outfile) fileName = TString(outfile);
    

    // create the task
    //===========================================================================
    AliAnalysisTaskDCArStudy *task = new AliAnalysisTaskDCArStudy(name);  
    if (!task) { return 0; }
    
    // configure the task
    //===========================================================================
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);    
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("default"));
    task->SetESDtrackCuts(0,AlidNdPtTools::CreateESDtrackCuts("TPCgeoNoDCAr"));    
    task->SetESDtrackCuts(1,AlidNdPtTools::CreateESDtrackCuts("TPCITSgeoNoDCAr"));    
        
    // attach the task to the manager and configure in and ouput
    //===========================================================================
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
  return task;
}