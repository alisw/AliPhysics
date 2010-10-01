//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Task for analysis
//  - reconstruction and MC output
// implementation file
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "TChain.h"
#include "TList.h"
#include "TH2F.h"

#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDtrackCuts.h"

#include "AliAnalysisTaskTotEt.h"
#include "AliAnalysisEtReconstructedPhos.h"
#include "AliAnalysisEtReconstructedEmcal.h"
#include "AliAnalysisEtMonteCarloPhos.h"
#include "AliAnalysisEtMonteCarloEmcal.h"

#include <iostream>

using namespace std;

ClassImp(AliAnalysisTaskTotEt)

//________________________________________________________________________
AliAnalysisTaskTotEt::AliAnalysisTaskTotEt(const char *name) :
        AliAnalysisTaskSE(name)
        ,fOutputList(0)
        ,fRecAnalysis(0)
        ,fMCAnalysis(0)
        ,fHistEtRecvsEtMC(0)
        ,fEsdtrackCutsTPC(0)
{
    // Constructor

    // select if we should use EMCal or PHOS class
    // PHOS by default, EMCal if name string contains EMC
    TString t(name);
    t.ToUpper();
    if (t.Contains("EMC")) {
      fRecAnalysis = new AliAnalysisEtReconstructedEmcal(); 
      fMCAnalysis = new AliAnalysisEtMonteCarloEmcal();
    }
    else {
      fRecAnalysis = new AliAnalysisEtReconstructedPhos(); 
      fMCAnalysis = new AliAnalysisEtMonteCarloPhos();
    }

    fRecAnalysis->Init();
    fMCAnalysis->Init();

    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #1 writes into a TH1 container

    DefineOutput(1, TList::Class());

}
AliAnalysisTaskTotEt::~AliAnalysisTaskTotEt(){//Destructor
  fOutputList->Clear();
  delete fOutputList;
  delete fRecAnalysis;
  delete fMCAnalysis;
  delete fEsdtrackCutsTPC;
}

//________________________________________________________________________
void AliAnalysisTaskTotEt::UserCreateOutputObjects()
{
    // Create histograms
    // Called once
    fMCAnalysis->CreateHistograms();
    fRecAnalysis->CreateHistograms();
    fOutputList = new TList;
    fOutputList->SetOwner();
    fRecAnalysis->FillOutputList(fOutputList);
    fMCAnalysis->FillOutputList(fOutputList);
    fHistEtRecvsEtMC = new TH2F("fHistEtRecvsEtMC", "Reconstructed E_{t} vs MC E_{t}", 1000, 0.000, 100, 1000, 0.0001, 100);
    fOutputList->Add(fHistEtRecvsEtMC);


    fEsdtrackCutsTPC = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fEsdtrackCutsTPC->SetName("fEsdTrackCutsTPCOnly");
    fOutputList->Add(fEsdtrackCutsTPC);
    if(fEsdtrackCutsTPC){
      fRecAnalysis->SetTPCOnlyTrackCuts( GetTPCOnlyTrackCuts());
      fMCAnalysis->SetTPCOnlyTrackCuts( GetTPCOnlyTrackCuts());
    }
    else{
      Printf("Error: no track cuts!");
    }

}

//________________________________________________________________________
void AliAnalysisTaskTotEt::UserExec(Option_t *)
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
        fMCAnalysis->AnalyseEvent(mcEvent);
    }

    fHistEtRecvsEtMC->Fill(fRecAnalysis->GetTotEtAcc(), fMCAnalysis->GetTotEt());

// Post output data.
    PostData(1, fOutputList);

}

//________________________________________________________________________
void AliAnalysisTaskTotEt::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query

    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) {
        printf("ERROR: Output list not available\n");
        return;
    }
}



