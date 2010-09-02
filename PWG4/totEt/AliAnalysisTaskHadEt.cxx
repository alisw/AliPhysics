//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TObjArray.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "TDatabasePDG.h"
#include "AliAnalysisTaskHadEt.h"

#include <iostream>
#include "AliStack.h"

using namespace std;

ClassImp(AliAnalysisTaskHadEt)



//________________________________________________________________________
AliAnalysisTaskHadEt::AliAnalysisTaskHadEt(const char *name) :
        AliAnalysisTaskSE(name)
        ,fESD(0)
        ,fOutputList(0)
        ,fRecAnalysis(0)
        ,fMCAnalysis(0)
        ,fHistEtRecvsEtMC(0)
        ,fTriggerSelection(false)
        ,fCount(0)
        ,fkPhotonPdg(22)
        ,fkProtonMass(.938)
        ,fPdgDB(0)
        ,fRecEventVars(0)
        ,fSimEventVars(0)
        ,esdtrackCutsITSTPC(0)
        ,esdtrackCutsTPC(0)
        ,esdtrackCutsITS(0)
{
    // Constructor

    fRecAnalysis = new AliAnalysisHadEtReconstructed();
    fRecAnalysis->Init();
    fMCAnalysis = new AliAnalysisHadEtMonteCarlo();
    fMCAnalysis->Init();
    if(!fPdgDB) fPdgDB = new TDatabasePDG();

    fTriggerSelection = false;
    fCount = 0;

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
    esdtrackCutsITSTPC = AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(selectPrimaries);
    esdtrackCutsITSTPC->SetName("fEsdTrackCuts");
    esdtrackCutsTPC = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    esdtrackCutsTPC->SetName("fEsdTrackCutsTPCOnly");
    //ITS stand alone cuts - similar to 2009 cuts but with only ITS hits required
    esdtrackCutsITS =  new AliESDtrackCuts;
    esdtrackCutsITS->SetName("fEsdTrackCutsITS");
    esdtrackCutsITS->SetRequireITSRefit(kTRUE);
    esdtrackCutsITS->SetRequireITSStandAlone(kTRUE);
    esdtrackCutsITS->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    esdtrackCutsITS->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");
    esdtrackCutsITS->SetMaxDCAToVertexZ(1.e6);
    esdtrackCutsITS->SetDCAToVertex2D(kFALSE);
    esdtrackCutsITS->SetRequireSigmaToVertex(kFALSE);
    esdtrackCutsITS->SetAcceptKinkDaughters(kFALSE);

    fOutputList->Add(esdtrackCutsITSTPC);
    fOutputList->Add(esdtrackCutsTPC);
    fOutputList->Add(esdtrackCutsITS);
    if(esdtrackCutsITSTPC && esdtrackCutsTPC){
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
{
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


