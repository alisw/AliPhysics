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

#include "AliProtonQAAnalysis.h"
#include "AliAnalysisTaskProtonsQA.h"

// Analysis task used for the QA of the (anti)proton analysis
// Author: Panos Cristakoglou

ClassImp(AliAnalysisTaskProtonsQA)

//________________________________________________________________________ 
AliAnalysisTaskProtonsQA::AliAnalysisTaskProtonsQA()
  : AliAnalysisTask(), fESD(0), fMC(0),
    fList0(0), fList1(0), fList2(0), fList3(0), fList4(0), fList5(0),
    fAnalysis(0) {
    //Dummy constructor
}

//________________________________________________________________________
AliAnalysisTaskProtonsQA::AliAnalysisTaskProtonsQA(const char *name) 
: AliAnalysisTask(name, ""), fESD(0), fMC(0),
  fList0(0), fList1(0), fList2(0), fList3(0), fList4(0), fList5(0),
  fAnalysis(0) {
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskProtonsQA::ConnectInputData(Option_t *) {
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
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
  fAnalysis = new AliProtonQAAnalysis();
  fAnalysis->SetRunMCAnalysis();
  fAnalysis->SetRunEfficiencyAnalysis(kTRUE,kTRUE); //kTRUE,kTRUE for eta-pT efficiencies and if the cuts should be used in the reco and pid efficiencies
  //fAnalysis->SetMCProcessId(13);//4: weak decay - 13: hadronic interaction
  //fAnalysis->SetMotherParticlePDGCode(3122);//3122: Lambda

  //Use of TPConly tracks
  /*fAnalysis->SetQAYPtBins(10, -0.5, 0.5, 12, 0.5, 0.9); //TPC only
  fAnalysis->UseTPCOnly();
  fAnalysis->SetTPCpid();
  fAnalysis->SetMinTPCClusters(100);
  fAnalysis->SetMaxChi2PerTPCCluster(2.2);
  fAnalysis->SetMaxCov11(0.5);
  fAnalysis->SetMaxCov22(0.5);
  fAnalysis->SetMaxCov33(0.5);
  fAnalysis->SetMaxCov44(0.5);
  fAnalysis->SetMaxCov55(0.5);
  //fAnalysis->SetMaxSigmaToVertexTPC(2.0);
  fAnalysis->SetMaxDCAXYTPC(1.5);
  fAnalysis->SetMaxDCAZTPC(1.5);*/

  //Use of HybridTPC tracks
  fAnalysis->SetQAYPtBins(20, -1.0, 1.0, 26, 0.2, 1.5); //HybridTPC
  fAnalysis->UseHybridTPC();
  fAnalysis->SetTPCpid();
  fAnalysis->SetMinTPCClusters(110);
  fAnalysis->SetMaxChi2PerTPCCluster(2.2);
  fAnalysis->SetMaxCov11(0.5);
  fAnalysis->SetMaxCov22(0.5);
  fAnalysis->SetMaxCov33(0.5);
  fAnalysis->SetMaxCov44(0.5);
  fAnalysis->SetMaxCov55(0.5);
  fAnalysis->SetMaxSigmaToVertex(2.0);
  //fAnalysis->SetMaxDCAXY(1.5);
  //fAnalysis->SetMaxDCAZ(1.5);
  //fAnalysis->SetPointOnITSLayer6();
  //fAnalysis->SetPointOnITSLayer5();
  //fAnalysis->SetPointOnITSLayer4();
  //fAnalysis->SetPointOnITSLayer3();
  fAnalysis->SetPointOnITSLayer2();
  fAnalysis->SetPointOnITSLayer1();
  fAnalysis->SetMinITSClusters(5);

  //Combined tracking
  /*fAnalysis->SetQAYPtBins(10, -0.5, 0.5, 12, 0.5, 0.9); //HybridTPC
  //fAnalysis->SetQAYPtBins(20, -1.0, 1.0, 27, 0.4, 3.1); //combined tracking
  fAnalysis->SetMinTPCClusters(110);
  fAnalysis->SetMaxChi2PerTPCCluster(2.2);
  fAnalysis->SetMaxCov11(0.5);
  fAnalysis->SetMaxCov22(0.5);
  fAnalysis->SetMaxCov33(0.5);
  fAnalysis->SetMaxCov44(0.5);
  fAnalysis->SetMaxCov55(0.5);
  fAnalysis->SetMaxSigmaToVertex(2.0);
  //fAnalysis->SetMaxDCAXY(2.0);
  //fAnalysis->SetMaxDCAZ(2.0);
  fAnalysis->SetTPCRefit();
  fAnalysis->SetPointOnITSLayer1();
  fAnalysis->SetPointOnITSLayer2();
  //fAnalysis->SetPointOnITSLayer3();
  //fAnalysis->SetPointOnITSLayer4();
  //fAnalysis->SetPointOnITSLayer5();
  //fAnalysis->SetPointOnITSLayer6();
  fAnalysis->SetMinITSClusters(1);
  fAnalysis->SetITSRefit();
  fAnalysis->SetESDpid();*/

  fAnalysis->SetPriorProbabilities(partFrac);

  fList0 = new TList();
  fList0 = fAnalysis->GetGlobalQAList();

  fList1 = new TList();
  fList1 = fAnalysis->GetPDGList();

  fList2 = new TList();
  fList2 = fAnalysis->GetMCProcessesList();

  fList3 = new TList();
  fList3 = fAnalysis->GetAcceptedCutList();

  fList4 = new TList();
  fList4 = fAnalysis->GetAcceptedDCAList();

  fList5 = new TList();
  fList5 = fAnalysis->GetEfficiencyQAList();
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
  
  fAnalysis->RunQAAnalysis(stack, fESD);
  fAnalysis->RunMCAnalysis(stack);
  fAnalysis->RunEfficiencyAnalysis(stack, fESD);

  // Post output data.
  PostData(0, fList0);
  PostData(1, fList1);
  PostData(2, fList2);
  PostData(3, fList3);
  PostData(4, fList4);
  PostData(5, fList5);
}      

//________________________________________________________________________
void AliAnalysisTaskProtonsQA::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
  
  fList0 = dynamic_cast<TList*> (GetOutputData(0));
  if (!fList0) {
    Printf("ERROR: fList not available");
    return;
  }
}


