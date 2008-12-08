#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TList.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliCFContainer.h"

#include "AliProtonAnalysis.h"
#include "AliAnalysisTaskProtons.h"

// Analysis task creating a the 2d y-p_t spectrum of p and antip
// Author: Panos Cristakoglou

ClassImp(AliAnalysisTaskProtons)

//________________________________________________________________________ 
AliAnalysisTaskProtons::AliAnalysisTaskProtons()
  : AliAnalysisTask(), fESD(0), fAOD(0), fMC(0), fAnalysisType("ESD"),
    fList(0), fAnalysis(0),
    fElectronFunction(0), fMuonFunction(0),
    fPionFunction(0), fKaonFunction(0), fProtonFunction(0),
    fFunctionUsed(kFALSE) {
  //Dummy constructor
                                                                                                   
}

//________________________________________________________________________
AliAnalysisTaskProtons::AliAnalysisTaskProtons(const char *name) 
: AliAnalysisTask(name, ""), fESD(0), fAOD(0), fMC(0), fAnalysisType("ESD"),
  fList(0), fAnalysis(0), 
  fElectronFunction(0), fMuonFunction(0),
  fPionFunction(0), fKaonFunction(0), fProtonFunction(0),
  fFunctionUsed(kFALSE) { 
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
    if(fAnalysisType == "ESD") {
      tree->SetBranchStatus("*", kFALSE);
      tree->SetBranchStatus("Tracks.*", kTRUE);
      
      AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
     
      if (!esdH) {
	Printf("ERROR: Could not get ESDInputHandler");
      } else
	fESD = esdH->GetEvent();
    }
    else if(fAnalysisType == "AOD") {
     AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
     
     if (!aodH) {
       Printf("ERROR: Could not get AODInputHandler");
     } else
       fAOD = aodH->GetEvent();
    }
    else if(fAnalysisType == "MC") {
     AliMCEventHandler* mcH = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
     if (!mcH) {
       Printf("ERROR: Could not retrieve MC event handler");
     }
     else
       fMC = mcH->MCEvent();
    }
    else
      Printf("Wrong analysis type: Only ESD, AOD and MC types are allowed!");
  }
}

//________________________________________________________________________
void AliAnalysisTaskProtons::CreateOutputObjects() {
  // Create histograms
  // Called once
  //Prior probabilities
  Double_t partFrac[5] = {0.01, 0.01, 0.85, 0.10, 0.05};
  
  //proton analysis object
  fAnalysis = new AliProtonAnalysis();
  fAnalysis->InitAnalysisHistograms(10,-1.0,1.0,26,0.5,3.1);

  if(fAnalysisType == "ESD") {
    //Use of TPConly tracks
    fAnalysis->UseTPCOnly();

    //TPC related cuts       
    fAnalysis->SetMinTPCClusters(50);
    fAnalysis->SetMaxChi2PerTPCCluster(3.5);
    fAnalysis->SetMaxCov11(2.0);
    fAnalysis->SetMaxCov22(2.0);
    fAnalysis->SetMaxCov33(0.5);
    fAnalysis->SetMaxCov44(0.5);
    fAnalysis->SetMaxCov55(2.0);
    fAnalysis->SetMaxSigmaToVertex(2.5);
    fAnalysis->SetTPCRefit();

    //ITS related cuts - to be used for the analysis of global tracking 
    //fAnalysis->SetMinITSClusters(5);
    //fAnalysis->SetITSRefit();
  }

  if(fFunctionUsed)
    fAnalysis->SetPriorProbabilityFunctions(fElectronFunction,
					    fMuonFunction,
					    fPionFunction,
					    fKaonFunction,
					    fProtonFunction);
  else
    fAnalysis->SetPriorProbabilities(partFrac);

  fList = new TList();
  fList->Add(fAnalysis->GetProtonYPtHistogram());
  fList->Add(fAnalysis->GetAntiProtonYPtHistogram());
  fList->Add(fAnalysis->GetEventHistogram());
  fList->Add(fAnalysis->GetProtonContainer());
  fList->Add(fAnalysis->GetAntiProtonContainer());
}

//________________________________________________________________________
void AliAnalysisTaskProtons::Exec(Option_t *) {
  // Main loop
  // Called for each event
  
  if(fAnalysisType == "ESD") {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }

    Printf("Proton ESD analysis task: There are %d tracks in this event", fESD->GetNumberOfTracks());
    fAnalysis->Analyze(fESD);
  }//ESD analysis              
  
  else if(fAnalysisType == "AOD") {
    if (!fAOD) {
      Printf("ERROR: fAOD not available");
      return;
    }
    
    Printf("Proton AOD analysis task: There are %d tracks in this event", fAOD->GetNumberOfTracks());
    fAnalysis->Analyze(fAOD);
  }//AOD analysis                                                                                                                                                                
  
  else if(fAnalysisType == "MC") {
    if (!fMC) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }

    AliStack* stack = fMC->Stack();
    if (!stack) {
      Printf("ERROR: Could not retrieve the stack");
      return;
    }
    Printf("Proton MC analysis task: There are %d primaries in this event", stack->GetNprimary());
    fAnalysis->Analyze(stack);
  }//MC analysis                      

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
    
  TCanvas *c1 = new TCanvas("c1","p-\bar{p}",200,0,800,400);
  c1->SetFillColor(10); c1->SetHighLightColor(10); c1->Divide(2,1);

  c1->cd(1)->SetLeftMargin(0.15); c1->cd(1)->SetBottomMargin(0.15);  
  if (fHistYPtProtons) fHistYPtProtons->DrawCopy("colz");
  c1->cd(2)->SetLeftMargin(0.15); c1->cd(2)->SetBottomMargin(0.15);  
  if (fHistYPtAntiProtons) fHistYPtAntiProtons->DrawCopy("colz");
}


