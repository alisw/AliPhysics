#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliBalance.h"

#include "AliAnalysisTaskBF.h"

// Analysis task for the BF code
// Authors: Panos Cristakoglou@cern.ch

ClassImp(AliAnalysisTaskBF)

//________________________________________________________________________
AliAnalysisTaskBF::AliAnalysisTaskBF(const char *name) 
  : AliAnalysisTask(name, ""), fESD(0), fAOD(0), fMC(0), fBalance(0) {
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
  TString gAnalysisLevel = fBalance->GetAnalysisLevel();

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } 
  else {
    if(gAnalysisLevel == "ESD") {
      AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
      
      if (!esdH) {
	Printf("ERROR: Could not get ESDInputHandler");
      } 
      else
	fESD = esdH->GetEvent();
    }//ESD
    else if(gAnalysisLevel == "AOD") {
      AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
      
      if (!aodH) {
	Printf("ERROR: Could not get AODInputHandler");
      } else
	fAOD = aodH->GetEvent();
    }//AOD
    else if(gAnalysisLevel == "MC") {
      AliMCEventHandler* mcH = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      if (!mcH) {
	Printf("ERROR: Could not retrieve MC event handler");
      }
      else
	fMC = mcH->MCEvent();
    }//MC
    else
      Printf("Wrong analysis type: Only ESD, AOD and MC types are allowed!");
  }
}

//________________________________________________________________________
void AliAnalysisTaskBF::CreateOutputObjects() {
  // Create histograms
  // Called once
  
}

//________________________________________________________________________
void AliAnalysisTaskBF::Exec(Option_t *) {
  // Main loop
  // Called for each event
  TString gAnalysisLevel = fBalance->GetAnalysisLevel();

  TObjArray *array = new TObjArray();
  
  //ESD analysis
  if(gAnalysisLevel == "ESD") {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }
    
    Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
    for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
      AliESDtrack* track = fESD->GetTrack(iTracks);
      if (!track) {
	Printf("ERROR: Could not receive track %d", iTracks);
	continue;
      }
      array->Add(track);
    } //track loop
  }//ESD analysis
  //AOD analysis
  else if(gAnalysisLevel == "AOD") {
    if (!fAOD) {
      Printf("ERROR: fAOD not available");
      return;
    }
    
    Printf("There are %d tracks in this event", fAOD->GetNumberOfTracks());
    for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
      AliAODtrack* track = fAOD->GetTrack(iTracks);
      if (!track) {
	Printf("ERROR: Could not receive track %d", iTracks);
	continue;
      }
      array->Add(track);
    } //track loop
  }//AOD analysis
  //MC analysis
  else if(gAnalysisLevel == "MC") {
    if (!fMC) {
      Printf("ERROR: fMC not available");
      return;
    }
    
    Printf("There are %d tracks in this event", fMC->GetNumberOfPrimaries());
    for (Int_t iTracks = 0; iTracks < fMC->GetNumberOfPrimaries(); iTracks++) {
      AliMCParticle* track = fMC->GetTrack(iTracks);
      if (!track) {
	Printf("ERROR: Could not receive particle %d", iTracks);
	continue;
      }
      array->Add(track);
    } //track loop
  }//MC analysis

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
