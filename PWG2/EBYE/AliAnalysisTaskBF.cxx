#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliBalance.h"

#include "AliAnalysisTaskBF.h"

// Analysis task for the BF code
// Authors: Panos Cristakoglou@cern.ch

ClassImp(AliAnalysisTaskBF)

//________________________________________________________________________
AliAnalysisTaskBF::AliAnalysisTaskBF(const char *name) 
  : AliAnalysisTask(name, ""), fESD(0), fBalance(0) {
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, AliBalance::Class());
}

//________________________________________________________________________
void AliAnalysisTaskBF::ConnectInputData(Option_t *) {
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } 
  else {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } 
    else
      fESD = esdH->GetEvent();
  }
}

//________________________________________________________________________
void AliAnalysisTaskBF::CreateOutputObjects() {
  // Create histograms
  // Called once

  fBalance = new AliBalance();
  fBalance->SetAnalysisType(1);
  fBalance->SetNumberOfBins(18);
  fBalance->SetInterval(-0.9,0.9);
  //cout<<fBalance->GetAnalysisType()<<endl;
}

//________________________________________________________________________
void AliAnalysisTaskBF::Exec(Option_t *) {
  // Main loop
  // Called for each event
  
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());

  TObjArray *array = new TObjArray();

  // Track loop to fill a pT spectrum
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    array->Add(track);
  } //track loop 
  fBalance->CalculateBalance(array);
  
  delete array;
  // Post output data.
  PostData(0, fBalance);
}      

//________________________________________________________________________
void AliAnalysisTaskBF::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  fBalance = dynamic_cast<AliBalance*> (GetOutputData(0));
  if (!fBalance) {
    Printf("ERROR: fBalance not available");
    return;
  }
  
  TGraphErrors *gr = fBalance->DrawBalance();
  gr->SetMarkerStyle(20);
  gr->Draw("AP");

  fBalance->PrintResults();
}
