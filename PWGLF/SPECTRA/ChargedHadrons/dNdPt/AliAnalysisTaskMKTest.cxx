#include "THnSparse.h"
#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliMultSelection.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisUtils.h"
#include "AliESDUtils.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AlidNdPtTools.h"
#include "AliAnalysisTaskMKTest.h"

class AliAnalysisTaskMKTest;    

using namespace std;

ClassImp(AliAnalysisTaskMKTest)

//_____________________________________________________________________________

AliAnalysisTaskMKTest::AliAnalysisTaskMKTest() 
    : AliAnalysisTaskSE()   
    , fEvent(0)
    , fESD(0)
    , fMCEvent(0)
    , fMCStack(0)
    , fMCHeader(0)
    , fMCGenHeader(0)
    , fESDtrackCuts0(0)
    , fOutputList(0)
    , fLogHist(0)
    , fHistDCA(0)    
{
    //constructor
}

//_____________________________________________________________________________

AliAnalysisTaskMKTest::AliAnalysisTaskMKTest(const char* name)
    : AliAnalysisTaskSE(name)
    , fEvent(0)
    , fESD(0)
    , fMCEvent(0)
    , fMCStack(0)
    , fMCHeader(0)
    , fMCGenHeader(0)
    , fESDtrackCuts0(0)
    , fOutputList(0)
    , fLogHist(0)
    , fHistDCA(0)    
{
    // constructor
     DefineInput(0, TChain::Class());    // define input
     DefineOutput(1, TList::Class());    // define ouptut 
}

//_____________________________________________________________________________

AliAnalysisTaskMKTest::~AliAnalysisTaskMKTest()
{
    // destructor
    if(fOutputList) { delete fOutputList; }
}

//_____________________________________________________________________________

void AliAnalysisTaskMKTest::UserCreateOutputObjects()
{
    // create output list
    fOutputList = new TList();         
    fOutputList->SetOwner(kTRUE);
    
    fLogHist = AlidNdPtTools::CreateLogHist("fLogHist");
    fOutputList->Add(fLogHist);

    AlidNdPtTools::AddAxis("pt");
    AlidNdPtTools::AddAxis("DCAxy",2000,-1,1);
    fHistDCA = AlidNdPtTools::CreateHist("fHistDCA");
    fOutputList->Add(fHistDCA);
    
    PostData(1, fOutputList);           // postdata 
   
}

//_____________________________________________________________________________
void AliAnalysisTaskMKTest::UserExec(Option_t *)
{           
    
    fEvent = InputEvent();
    if (!fEvent) { Log("noEvent"); return; }
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());   
    if (!fESD) { Log("noESD"); return; }

    fMCEvent = MCEvent();
    if (fMCEvent) { Log("MC"); }
      fMCHeader = fMCEvent->Header();
      if (!fMCHeader) { Log("noMCHeader"); }
      fMCGenHeader = fMCHeader->GenEventHeader();
      if (!fMCGenHeader) { 
          Log("noMCGenHeader");           
      } else { 
          TString s = "mcHeader=";
          s += fMCGenHeader->GetName();
          Log(s.Data());
      }
      fMCStack = fMCEvent->Stack();    
      if (!fMCStack) { Log("noMCStack"); }    
      
    // first loop for nacc and maxpt
    Int_t nTracks = fESD->GetNumberOfTracks();
    AliESDtrack* track = 0;    
    Double_t maxpt = 0; 
    Int_t nacc = 0;

    for(Int_t iTrack = 0; iTrack < nTracks; iTrack++) {     // track loop
        track = static_cast<AliESDtrack*>(fESD->GetTrack(iTrack));
        if (!track) {
            Log("noESDtrackLoop0");
            continue;
        }    
        if (! fESDtrackCuts0 ) { continue; } 
        if (! fESDtrackCuts0->AcceptTrack(track) ) { continue; } 
        nacc++;
        if (track->Pt() > maxpt) { maxpt = track->Pt(); }
    }

    // second loop for track hist filling
    for(Int_t iTrack = 0; iTrack < nTracks; iTrack++) {     // track loop
        track = static_cast<AliESDtrack*>(fESD->GetTrack(iTrack));
        if (!track) {
            Log("noESDtrack");
            continue;
        }
        
          Float_t b[2];
          Float_t bCov[3];
       
          track->GetImpactParameters(b,bCov);
        
        Float_t dcar = b[0];
        Float_t dcaz = b[1];

        
        Double_t pt = track->Pt();
        Double_t eta = track->Eta();
        Double_t q   = track->GetSign();        

        if (! fESDtrackCuts0 ) { continue; } 
        if (! fESDtrackCuts0->AcceptTrack(track) ) { continue; } 
        AlidNdPtTools::FillHist(fHistDCA, pt, dcar);
        }
      
} 


//_____________________________________________________________________________

void AliAnalysisTaskMKTest::FinishTaskOutput()
{
    // finish task output
    PostData(1, fOutputList); 
}

//_____________________________________________________________________________

void AliAnalysisTaskMKTest::Terminate(Option_t *)
{
    // terminate
}

//____________________________________________________________________________

void AliAnalysisTaskMKTest::Log(const char* name) 
{ 
    AlidNdPtTools::Log(fLogHist, name);   
}

//____________________________________________________________________________