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
#include "AliAnalysisTaskTrackCuts.h"

class AliAnalysisTaskTrackCuts;

namespace {
using namespace std;
}
/// \cond CLASSIMP
ClassImp(AliAnalysisTaskTrackCuts)
/// \endcond
//_____________________________________________________________________________

AliAnalysisTaskTrackCuts::AliAnalysisTaskTrackCuts()
    : AliAnalysisTaskMKBase()
{
    // default contructor
}

//_____________________________________________________________________________

AliAnalysisTaskTrackCuts::AliAnalysisTaskTrackCuts(const char* name)
    : AliAnalysisTaskMKBase(name)
{
    // constructor
}

//_____________________________________________________________________________

void AliAnalysisTaskTrackCuts::Terminate(Option_t *)
{
    // terminate
    fESDtrackCutsM->SaveHistograms(0);
    cout<<"AliAnalysisTaskTrackCuts terminated."<<endl;
}

//_____________________________________________________________________________

AliAnalysisTaskTrackCuts* AliAnalysisTaskTrackCuts::AddTaskTrackCuts(const char* name, const char* outfile, int _cutMode)
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskTrackCuts", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskTrackCuts", "This task requires an input event handler");
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
    AliAnalysisTaskTrackCuts *task = new AliAnalysisTaskTrackCuts(name);
    if (!task) { return 0; }
    
    // configure the task
    //===========================================================================
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta08", _cutMode, true));
    // attach the task to the manager and configure in and ouput
    //===========================================================================
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
  return task;
}
