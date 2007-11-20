#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "PWG2spectra/SPECTRA/AliProtonAnalysis.h"
#include "AliAnalysisTaskProtons.h"

// Analysis task creating a the 2d y-p_t spectrum of p and antip
// Author: Panos Cristakoglou

ClassImp(AliAnalysisTaskProtons)

//________________________________________________________________________
AliAnalysisTaskProtons::AliAnalysisTaskProtons(const char *name) 
  : AliAnalysisTask(name, ""), fESD(0), fList(0), fAnalysis(0) { 
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskProtons::ConnectInputData(Option_t *) {
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    // Disable all branches and enable only the needed ones
    // The next two lines are different when data produced as AliESDEvent is read
    tree->SetBranchStatus("*", kFALSE);
    tree->SetBranchStatus("fTracks.*", kTRUE);

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESD = esdH->GetEvent();
  }
}

//________________________________________________________________________
void AliAnalysisTaskProtons::CreateOutputObjects() {
  // Create histograms
  // Called once
  Double_t partFrac[5] = {0.01, 0.01, 0.85, 0.10, 0.05};

  fAnalysis = new AliProtonAnalysis();
  fAnalysis->InitHistograms(10,-1.0,1.0,30,0.1,3.1);
  fAnalysis->SetMinTPCClusters(50);
  fAnalysis->SetMinITSClusters(1);
  fAnalysis->SetMaxChi2PerTPCCluster(3.5);
  fAnalysis->SetMaxCov11(2.0);
  fAnalysis->SetMaxCov22(2.0);
  fAnalysis->SetMaxCov33(0.5);
  fAnalysis->SetMaxCov44(0.5);
  fAnalysis->SetMaxCov55(2.0);
  fAnalysis->SetMaxSigmaToVertex(3.);
  fAnalysis->SetITSRefit();
  fAnalysis->SetTPCRefit();
  fAnalysis->SetPriorProbabilities(partFrac);
}

//________________________________________________________________________
void AliAnalysisTaskProtons::Exec(Option_t *) {
  // Main loop
  // Called for each event
  
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }
  
  Printf("Task2: There are %d tracks in this event", fESD->GetNumberOfTracks());
  fAnalysis->Analyze(fESD);
  fList = new TList();
  fList->Add(fAnalysis->GetProtonYPtHistogram()); 
  fList->Add(fAnalysis->GetAntiProtonYPtHistogram()); 

  // Post output data.
  PostData(0, fList);
}      

//________________________________________________________________________
void AliAnalysisTaskProtons::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
  
  fList = dynamic_cast<TList*> (GetOutputData(0));
  if (!fList) {
    Printf("ERROR: fList not available");
    return;
  }
   
  TH2F *fHistYPtProtons = (TH2F *)fList->At(0);
  TH2F *fHistYPtAntiProtons = (TH2F *)fList->At(1);
    
  TCanvas *c2 = new TCanvas("c2","p-\bar{p}",200,0,800,400);
  c2->SetFillColor(10); c2->SetHighLightColor(10); c2->Divide(2,1);

  c2->cd(1)->SetLeftMargin(0.15); c2->cd(1)->SetBottomMargin(0.15);  
  if (fHistYPtProtons) fHistYPtProtons->DrawCopy("colz");
  c2->cd(2)->SetLeftMargin(0.15); c2->cd(2)->SetBottomMargin(0.15);  
  if (fHistYPtAntiProtons) fHistYPtAntiProtons->DrawCopy("colz");
}


