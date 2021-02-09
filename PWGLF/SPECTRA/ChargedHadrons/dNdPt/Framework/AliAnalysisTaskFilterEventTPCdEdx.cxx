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
#include "AliAnalysisTaskFilterEventTPCdEdx.h"

class AliAnalysisTaskFilterEventTPCdEdx;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskFilterEventTPCdEdx)
/// \endcond
//_____________________________________________________________________________

AliAnalysisTaskFilterEventTPCdEdx::AliAnalysisTaskFilterEventTPCdEdx() 
    : AliAnalysisTaskMKBase()
    , fHistdEdx(0)
    , fesdTreeFiltered(0)
{
    // default contructor
}

//_____________________________________________________________________________

AliAnalysisTaskFilterEventTPCdEdx::AliAnalysisTaskFilterEventTPCdEdx(const char* name) 
    : AliAnalysisTaskMKBase(name)
    , fHistdEdx(0)    
    , fesdTreeFiltered(0)
{
    // constructor
    DefineOutput(2, TTree::Class());
}

//_____________________________________________________________________________

AliAnalysisTaskFilterEventTPCdEdx::~AliAnalysisTaskFilterEventTPCdEdx()
{
    // destructor
}

//_____________________________________________________________________________

void AliAnalysisTaskFilterEventTPCdEdx::AddOutput()
{    
    // pt:eta:q:p:dedxsignal
    AddAxis("pT","pt");          
    AddAxis("eta","#eta", 20, -1, 1);
    AddAxis("Q",2,-1.5,1.5);
    //AddAxis("p","p (GeV/c)","pt");          
    AddAxis("ld(p)",2000,-1,1);       
    AddAxis("TPCsignal",2000,0,2000);        
    fHistdEdx = CreateHist("fHistdEdx");
    fOutputList->Add(fHistdEdx);    
    OpenFile(2);  

    //
    // Create tree
    fesdTreeFiltered = new TTree("esdTree", "Tree with filtered ESDs");  
  
    PostData(2, fesdTreeFiltered);
}

//_____________________________________________________________________________

Bool_t AliAnalysisTaskFilterEventTPCdEdx::IsEventSelected()
{
    return fIsAcceptedAliEventCuts;
}

//_____________________________________________________________________________

void AliAnalysisTaskFilterEventTPCdEdx::AnaEvent()
{   
   LoopOverAllTracks();
      
   AliESDEvent* esdIn = fESD;
   TTree* esdTreeOut = fesdTreeFiltered;
   Int_t iev = fEventNumberInFile;
   if (esdIn) return;
   

    esdTreeOut->Fill();


   PostData(2, fesdTreeFiltered);
   
}

//_____________________________________________________________________________

void AliAnalysisTaskFilterEventTPCdEdx::AnaTrack(Int_t flag)
{
    if (!fAcceptTrack[0]) return;
    // pt:eta:q:p:dedxsignal
    FillHist(fHistdEdx, fPt, fEta, fChargeSign, TMath::Log10(fP), fESDTrack->GetTPCsignal());
}

//_____________________________________________________________________________

AliAnalysisTaskFilterEventTPCdEdx* AliAnalysisTaskFilterEventTPCdEdx::AddTaskFilterEventTPCdEdx(const char* name, const char* treefile, const char* outfile) 
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskFilterEventTPCdEdx", "No analysis manager to connect to.");
        return 0;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskFilterEventTPCdEdx", "This task requires an input event handler");
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
    
    TString fileNameTree = AliAnalysisManager::GetCommonFileName();
    fileNameTree += ":";
    fileNameTree += name; 
    if (treefile) { // if a finename is given, use that one
        fileNameTree = TString(treefile);        
    }    
  
    // create the task
    //===========================================================================
    AliAnalysisTaskFilterEventTPCdEdx *task = new AliAnalysisTaskFilterEventTPCdEdx(name);  
    if (!task) { return 0; }
    
    // configure the task
    //===========================================================================
    task->SelectCollisionCandidates(AliVEvent::kAnyINT);    
    task->SetESDtrackCutsM(AlidNdPtTools::CreateESDtrackCuts("defaultEta08"));
    task->SetESDtrackCuts(0,AlidNdPtTools::CreateESDtrackCuts("TPConlyMinimalEta10"));
    task->GetESDtrackCuts(0)->SetMinNClustersTPC(70);    
    task->GetESDtrackCuts(0)->SetMaxDCAToVertexZ(20);
    task->GetESDtrackCuts(0)->SetMaxDCAToVertexXY(20);    
    task->GetESDtrackCuts(0)->SetRequireITSRefit(kTRUE);
    
    // attach the task to the manager and configure in and ouput
    //===========================================================================
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    
    // output container 1 for histograms    
    mgr->ConnectOutput(task, 1, mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
    // output container 2 for filtered esd event tree
    mgr->ConnectOutput(task, 2, mgr->CreateContainer("filteredESDs", TTree::Class(), AliAnalysisManager::kOutputContainer, fileNameTree.Data()));

    
  return task;
}