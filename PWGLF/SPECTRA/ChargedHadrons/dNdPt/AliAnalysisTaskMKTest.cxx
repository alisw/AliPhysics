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
#include "AliAnalysisTaskMKTest.h"

class AliAnalysisTaskMKTest;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskMKTest)
/// \endcond
//_____________________________________________________________________________

AliAnalysisTaskMKTest::AliAnalysisTaskMKTest() 
    : AliAnalysisTaskMKBase()
    , fHistPt(0)
{
    // default contructor
}

//_____________________________________________________________________________

AliAnalysisTaskMKTest::AliAnalysisTaskMKTest(const char* name) 
    : AliAnalysisTaskMKBase(name)
    , fHistPt(0)    
{
    // constructor
}

//_____________________________________________________________________________

AliAnalysisTaskMKTest::~AliAnalysisTaskMKTest()
{
    // destructor
}

//_____________________________________________________________________________

void AliAnalysisTaskMKTest::AddOutput()
{    
    AddAxis("pT","pt");    
    AddAxis("pT,inner","pt");    
    AddAxis("pT,TPCinner","pt");    
    AddAxis("pT,MC","pt");    
    fHistPt = CreateHist("fHistPt");
    fOutputList->Add(fHistPt);
    
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskMKTest::IsEventSelected()
{
    return fIsAcceptedAliEventCuts;
}

//_____________________________________________________________________________

void AliAnalysisTaskMKTest::AnaEvent()
{   
   LoopOverAllTracks();
}

//_____________________________________________________________________________

void AliAnalysisTaskMKTest::AnaTrack(Int_t flag)
{
    if (!fAcceptTrack[0]) return;
    FillHist(fHistPt, fPt, fPtInner, fPtInnerTPC, fMCPt);
}

//_____________________________________________________________________________

AliAnalysisTaskMKTest* AliAnalysisTaskMKTest::AddTaskMKTest(const char* name, const char* outfile) 
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskMKTest", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskMKTest", "This task requires an input event handler");
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
    AliAnalysisTaskMKTest *task = new AliAnalysisTaskMKTest(name);  
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