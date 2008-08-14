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
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "PWG2spectra/SPECTRA/AliProtonAnalysis.h"
#include "AliAnalysisTaskProtonsQA.h"

// Analysis task used for the QA of the (anti)proton analysis
// Author: Panos Cristakoglou

ClassImp(AliAnalysisTaskProtonsQA)

//________________________________________________________________________ 
AliAnalysisTaskProtonsQA::AliAnalysisTaskProtonsQA()
  : AliAnalysisTask(), fESD(0), fMC(0),
    fList(0), fAnalysis(0) {
    //Dummy constructor
                                                                                                   
}

//________________________________________________________________________
AliAnalysisTaskProtonsQA::AliAnalysisTaskProtonsQA(const char *name) 
: AliAnalysisTask(name, ""), fESD(0), fMC(0),
  fList(0), fAnalysis(0) {
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskProtonsQA::ConnectInputData(Option_t *) {
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    //tree->SetBranchStatus("*", kFALSE);
    //tree->SetBranchStatus("Tracks.*", kTRUE);
      
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
     
    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESD = esdH->GetEvent();
  }

  AliMCEventHandler* mcH = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!mcH) {
    Printf("ERROR: Could not retrieve MC event handler");
  }
  else
    fMC = mcH->MCEvent();
}

//________________________________________________________________________
void AliAnalysisTaskProtonsQA::CreateOutputObjects() {
  // Create histograms
  // Called once
  //Prior probabilities
  Double_t partFrac[5] = {0.01, 0.01, 0.85, 0.10, 0.05};
  
  //proton analysis object
  fAnalysis = new AliProtonAnalysis();
  fAnalysis->SetQAOn();
  fAnalysis->InitQA();

  //Use of TPConly tracks
  /*fAnalysis->SetQAYPtBins(10, -0.5, 0.5, 12, 0.5, 0.9); //TPC only
  fAnalysis->UseTPCOnly();
  fAnalysis->SetMinTPCClusters(50);
  fAnalysis->SetMaxChi2PerTPCCluster(3.5);
  fAnalysis->SetMaxCov11(2.0);
  fAnalysis->SetMaxCov22(2.0);
  fAnalysis->SetMaxCov33(0.5);
  fAnalysis->SetMaxCov44(0.5);
  fAnalysis->SetMaxCov55(2.0);
  fAnalysis->SetMaxSigmaToVertexTPC(2.5);
  fAnalysis->SetTPCRefit();
  fAnalysis->SetTPCpid();*/

  //Combined tracking
  fAnalysis->SetQAYPtBins(20, -1.0, 1.0, 27, 0.4, 3.1); //combined tracking
  fAnalysis->SetMinTPCClusters(50);
  fAnalysis->SetMaxChi2PerTPCCluster(3.5);
  fAnalysis->SetMaxCov11(2.0);
  fAnalysis->SetMaxCov22(2.0);
  fAnalysis->SetMaxCov33(0.5);
  fAnalysis->SetMaxCov44(0.5);
  fAnalysis->SetMaxCov55(2.0);
  fAnalysis->SetMaxSigmaToVertexTPC(2.5);
  fAnalysis->SetTPCRefit();
  //ITS related cuts - to be used in the case of the analysis of global tracks
  fAnalysis->SetMinITSClusters(5);
  fAnalysis->SetITSRefit();
  fAnalysis->SetESDpid();

  fAnalysis->SetPriorProbabilities(partFrac);

  fList = new TList();
  fList = fAnalysis->GetGlobalQAList();
}

//________________________________________________________________________
void AliAnalysisTaskProtonsQA::Exec(Option_t *) {
  // Main loop
  // Called for each event
  
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  if (!fMC) {
    Printf("ERROR: Could not retrieve MC event");
    return;
  }
  
  AliStack* stack = fMC->Stack();
  if (!stack) {
    Printf("ERROR: Could not retrieve the stack");
    return;
  }
  
  fAnalysis->RunQA(stack, fESD);
    
  // Post output data.
  PostData(0, fList);
}      

//________________________________________________________________________
void AliAnalysisTaskProtonsQA::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
  
  fList = dynamic_cast<TList*> (GetOutputData(0));
  if (!fList) {
    Printf("ERROR: fList not available");
    return;
  }
}


