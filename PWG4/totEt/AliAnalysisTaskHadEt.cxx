//_________________________________________________________________________
//  Utility Class for transverse energy studies; charged hadrons
//  Task for analysis
//  - reconstruction and MC output
// implementation file
//
//Created by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
//_________________________________________________________________________

#include "TChain.h"
#include "TList.h"
#include "TH2F.h"

#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDtrackCuts.h"

#include "AliAnalysisTaskHadEt.h"
#include "AliAnalysisHadEtReconstructed.h"
#include "AliAnalysisHadEtMonteCarlo.h"

#include <iostream>

using namespace std;

ClassImp(AliAnalysisTaskHadEt)



//________________________________________________________________________
AliAnalysisTaskHadEt::AliAnalysisTaskHadEt(const char *name) :
        AliAnalysisTaskSE(name)
        ,fOutputList(0)
        ,fRecAnalysis(0)
        ,fMCAnalysis(0)
        ,fHistEtRecvsEtMC(0)
        ,fEsdtrackCutsITSTPC(0)
        ,fEsdtrackCutsTPC(0)
        ,fEsdtrackCutsITS(0)
{
    // Constructor

    fRecAnalysis = new AliAnalysisHadEtReconstructed();
    fRecAnalysis->Init();
    fMCAnalysis = new AliAnalysisHadEtMonteCarlo();
    fMCAnalysis->Init();

    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #1 writes into a TH1 container

    DefineOutput(1, TList::Class());

}


//________________________________________________________________________
void AliAnalysisTaskHadEt::UserCreateOutputObjects()
{
    // Create histograms
    // Called once
  fOutputList = new TList;
  fMCAnalysis->SetHistoList(fOutputList);
  fRecAnalysis->SetHistoList(fOutputList);
  fMCAnalysis->CreateHistograms();
  fRecAnalysis->CreateHistograms();
  fRecAnalysis->FillOutputList();
  fMCAnalysis->FillOutputList();


    Bool_t selectPrimaries=kTRUE;
    fEsdtrackCutsITSTPC = AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(selectPrimaries);
    fEsdtrackCutsITSTPC->SetName("fEsdTrackCuts");
    fEsdtrackCutsTPC = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fEsdtrackCutsTPC->SetName("fEsdTrackCutsTPCOnly");
    //ITS stand alone cuts - similar to 2009 cuts but with only ITS hits required
    fEsdtrackCutsITS =  new AliESDtrackCuts;
    fEsdtrackCutsITS->SetName("fEsdTrackCutsITS");
    fEsdtrackCutsITS->SetRequireITSRefit(kTRUE);
    fEsdtrackCutsITS->SetRequireITSStandAlone(kTRUE);
    fEsdtrackCutsITS->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    fEsdtrackCutsITS->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");
    fEsdtrackCutsITS->SetMaxDCAToVertexZ(1.e6);
    fEsdtrackCutsITS->SetDCAToVertex2D(kFALSE);
    fEsdtrackCutsITS->SetRequireSigmaToVertex(kFALSE);
    fEsdtrackCutsITS->SetAcceptKinkDaughters(kFALSE);

    fOutputList->Add(fEsdtrackCutsITSTPC);
    fOutputList->Add(fEsdtrackCutsTPC);
    fOutputList->Add(fEsdtrackCutsITS);
    if(fEsdtrackCutsITSTPC && fEsdtrackCutsTPC){
      fRecAnalysis->SetITSTrackCuts( GetITSTrackCuts());
      fMCAnalysis->SetITSTrackCuts( GetITSTrackCuts());
      fRecAnalysis->SetTPCITSTrackCuts( GetTPCITSTrackCuts());
      fMCAnalysis->SetTPCITSTrackCuts( GetTPCITSTrackCuts());
      fRecAnalysis->SetTPCOnlyTrackCuts( GetTPCOnlyTrackCuts());
      fMCAnalysis->SetTPCOnlyTrackCuts( GetTPCOnlyTrackCuts());
      //add ITS stuff!
    }
    else{
      Printf("Error: no track cuts!");
    }

}

//________________________________________________________________________
void AliAnalysisTaskHadEt::UserExec(Option_t *)
{ // execute method
  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
if (!event) {
  Printf("ERROR: Could not retrieve event");
  return;
 }

fRecAnalysis->AnalyseEvent(event);

AliMCEvent* mcEvent = MCEvent();
if (mcEvent)
  {
    ((AliAnalysisHadEtMonteCarlo*)fMCAnalysis)->AnalyseEvent((AliVEvent*)mcEvent,(AliVEvent*)event);
  }

// Post output data.
PostData(1, fOutputList);
 
}

//________________________________________________________________________
void AliAnalysisTaskHadEt::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query

    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) {
        printf("ERROR: Output list not available\n");
        return;
    }
}


