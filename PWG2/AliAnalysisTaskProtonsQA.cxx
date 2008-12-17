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
#include "AliVertexerTracks.h"
#include "AliESDVertex.h"

#include "AliProtonQAAnalysis.h"
#include "AliAnalysisTaskProtonsQA.h"

// Analysis task used for the QA of the (anti)proton analysis
// Author: Panos Cristakoglou

ClassImp(AliAnalysisTaskProtonsQA)
  
//________________________________________________________________________ 
AliAnalysisTaskProtonsQA::AliAnalysisTaskProtonsQA()
  : AliAnalysisTask(), fESD(0), fMC(0),
    fList0(0), fList1(0), fList2(0), fList3(0), fList4(0), fList5(0),
    fProtonQAAnalysis(0), 
    fTriggerMode(kMB2), fProtonAnalysisMode(kTPC),
    fVxMax(0), fVyMax(0), fVzMax(0) {
  //Dummy constructor
}

//________________________________________________________________________
AliAnalysisTaskProtonsQA::AliAnalysisTaskProtonsQA(const char *name) 
  : AliAnalysisTask(name, ""), fESD(0), fMC(0),
    fList0(0), fList1(0), fList2(0), fList3(0), fList4(0), fList5(0),
    fProtonQAAnalysis(0),
    fTriggerMode(kMB2), fProtonAnalysisMode(kTPC),
    fVxMax(0), fVyMax(0), fVzMax(0) {
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
  fProtonQAAnalysis = new AliProtonQAAnalysis();
  fProtonQAAnalysis->SetRunMCAnalysis();
  fProtonQAAnalysis->SetRunEfficiencyAnalysis(kTRUE,kFALSE); //kTRUE,kTRUE for eta-pT efficiencies and if the cuts should be used in the reco and pid efficiencies
  //fProtonQAAnalysis->SetMCProcessId(13);//4: weak decay - 13: hadronic interaction
  //fProtonQAAnalysis->SetMotherParticlePDGCode(3122);//3122: Lambda

  //Use of TPConly tracks
  if(fProtonAnalysisMode == kTPC) {
    fProtonQAAnalysis->SetQAYPtBins(10, -0.5, 0.5, 12, 0.5, 0.9); //TPC only
    fProtonQAAnalysis->UseTPCOnly();
    //fProtonQAAnalysis->SetTPCpid();
    fProtonQAAnalysis->SetMinTPCClusters(100);
    fProtonQAAnalysis->SetMaxChi2PerTPCCluster(2.2);
    fProtonQAAnalysis->SetMaxCov11(0.5);
    fProtonQAAnalysis->SetMaxCov22(0.5);
    fProtonQAAnalysis->SetMaxCov33(0.5);
    fProtonQAAnalysis->SetMaxCov44(0.5);
    fProtonQAAnalysis->SetMaxCov55(0.5);
    fProtonQAAnalysis->SetMaxSigmaToVertexTPC(2.0);
    //fProtonQAAnalysis->SetMaxDCAXYTPC(1.5);
    //fProtonQAAnalysis->SetMaxDCAZTPC(1.5);
  }
  //Use of HybridTPC tracks
  else if(fProtonAnalysisMode == kHybrid) {
    fProtonQAAnalysis->SetQAYPtBins(10, -0.5, 0.5, 12, 0.5, 0.9); //HybridTPC
    fProtonQAAnalysis->UseHybridTPC();
    fProtonQAAnalysis->SetTPCpid();
    fProtonQAAnalysis->SetMinTPCClusters(110);
    fProtonQAAnalysis->SetMaxChi2PerTPCCluster(2.2);
    fProtonQAAnalysis->SetMaxCov11(0.5);
    fProtonQAAnalysis->SetMaxCov22(0.5);
    fProtonQAAnalysis->SetMaxCov33(0.5);
    fProtonQAAnalysis->SetMaxCov44(0.5);
    fProtonQAAnalysis->SetMaxCov55(0.5);
    fProtonQAAnalysis->SetMaxSigmaToVertex(2.0);
    /*fProtonQAAnalysis->SetMaxDCAXY(1.5);
      fProtonQAAnalysis->SetMaxDCAZ(1.5);*/
    fProtonQAAnalysis->SetPointOnITSLayer6();
    fProtonQAAnalysis->SetPointOnITSLayer5();
  //fProtonQAAnalysis->SetPointOnITSLayer4();
  //fProtonQAAnalysis->SetPointOnITSLayer3();
    fProtonQAAnalysis->SetPointOnITSLayer2();
    fProtonQAAnalysis->SetPointOnITSLayer1();
    fProtonQAAnalysis->SetMinITSClusters(5);
  }
  //Combined tracking
  else if(fProtonAnalysisMode == kGlobal) {
    fProtonQAAnalysis->SetQAYPtBins(20, -1.0, 1.0, 27, 0.4, 3.1); //combined tracking
    fProtonQAAnalysis->SetMinTPCClusters(110);
    fProtonQAAnalysis->SetMaxChi2PerTPCCluster(2.2);
    fProtonQAAnalysis->SetMaxCov11(0.5);
    fProtonQAAnalysis->SetMaxCov22(0.5);
    fProtonQAAnalysis->SetMaxCov33(0.5);
    fProtonQAAnalysis->SetMaxCov44(0.5);
    fProtonQAAnalysis->SetMaxCov55(0.5);
    fProtonQAAnalysis->SetMaxSigmaToVertex(2.0);
    //fProtonQAAnalysis->SetMaxDCAXY(2.0);
    //fProtonQAAnalysis->SetMaxDCAZ(2.0);
    fProtonQAAnalysis->SetTPCRefit();
    fProtonQAAnalysis->SetPointOnITSLayer1();
    fProtonQAAnalysis->SetPointOnITSLayer2();
    //fProtonQAAnalysis->SetPointOnITSLayer3();
    //fProtonQAAnalysis->SetPointOnITSLayer4();
    fProtonQAAnalysis->SetPointOnITSLayer5();
    fProtonQAAnalysis->SetPointOnITSLayer6();
    fProtonQAAnalysis->SetMinITSClusters(5);
    fProtonQAAnalysis->SetITSRefit();
    fProtonQAAnalysis->SetESDpid();
  }

  fProtonQAAnalysis->SetPriorProbabilities(partFrac);
  
  fList0 = new TList();
  fList0 = fProtonQAAnalysis->GetGlobalQAList();

  fList1 = new TList();
  fList1 = fProtonQAAnalysis->GetPDGList();

  fList2 = new TList();
  fList2 = fProtonQAAnalysis->GetMCProcessesList();

  fList3 = new TList();
  fList3 = fProtonQAAnalysis->GetAcceptedCutList();

  fList4 = new TList();
  fList4 = fProtonQAAnalysis->GetAcceptedDCAList();

  fList5 = new TList();
  fList5 = fProtonQAAnalysis->GetEfficiencyQAList();
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
  
  if(IsEventTriggered(fESD,fTriggerMode)) {
    const AliESDVertex *vertex = GetVertex(fESD,fProtonAnalysisMode,
					   fVxMax,fVyMax,fVzMax);
    if(vertex) {
      fProtonQAAnalysis->RunQAAnalysis(stack, fESD);
      fProtonQAAnalysis->RunMCAnalysis(stack);
      fProtonQAAnalysis->RunEfficiencyAnalysis(stack, fESD);
    }//accepted vertex
  }//triggered event
  
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
    Printf("ERROR: fList0 not available");
    return;
  }
  fList1 = dynamic_cast<TList*> (GetOutputData(1));
  if (!fList1) {
    Printf("ERROR: fList1 not available");
    return;
  }
  fList2 = dynamic_cast<TList*> (GetOutputData(2));
  if (!fList2) {
    Printf("ERROR: fList2 not available");
    return;
  }
  fList3 = dynamic_cast<TList*> (GetOutputData(3));
  if (!fList3) {
    Printf("ERROR: fList3 not available");
    return;
  }
  fList4 = dynamic_cast<TList*> (GetOutputData(4));
  if (!fList4) {
    Printf("ERROR: fList4 not available");
    return;
  }
  fList5 = dynamic_cast<TList*> (GetOutputData(5));
  if (!fList5) {
    Printf("ERROR: fList5 not available");
    return;
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskProtonsQA::IsEventTriggered(const AliESDEvent *esd, 
						  TriggerMode trigger) {
  // check if the event was triggered
  ULong64_t triggerMask = esd->GetTriggerMask();
  
  // definitions from p-p.cfg
  ULong64_t spdFO = (1 << 14);
  ULong64_t v0left = (1 << 11);
  ULong64_t v0right = (1 << 12);

  switch (trigger) {
  case kMB1: {
    if (triggerMask & spdFO || ((triggerMask & v0left) || (triggerMask & v0right)))
      return kTRUE;
    break;
  }
  case kMB2: {
    if (triggerMask & spdFO && ((triggerMask & v0left) || (triggerMask & v0right)))
      return kTRUE;
    break;
  }
  case kSPDFASTOR: {
    if (triggerMask & spdFO)
      return kTRUE;
    break;
  }
  }//switch

  return kFALSE;
}

//________________________________________________________________________
const AliESDVertex* AliAnalysisTaskProtonsQA::GetVertex(AliESDEvent* esd, 
							AnalysisMode mode,
							Double_t gVxMax,
							Double_t gVyMax,
							Double_t gVzMax) {
  // Get the vertex from the ESD and returns it if the vertex is valid
  // Second argument decides which vertex is used (this selects
  // also the quality criteria that are applied)
  const AliESDVertex* vertex = 0;
  if(mode == kHybrid)
    vertex = esd->GetPrimaryVertexSPD();
  else if(mode == kTPC){
    Double_t kBz = esd->GetMagneticField();
    AliVertexerTracks vertexer(kBz);
    vertexer.SetTPCMode();
    AliESDVertex *vTPC = vertexer.FindPrimaryVertex(esd);
    esd->SetPrimaryVertexTPC(vTPC);
    for (Int_t i=0; i<esd->GetNumberOfTracks(); i++) {
      AliESDtrack *t = esd->GetTrack(i);
      t->RelateToVertexTPC(vTPC, kBz, kVeryBig);
    }
    delete vTPC;
    vertex = esd->GetPrimaryVertexTPC();
  }
  else if(mode == kGlobal)
    vertex = esd->GetPrimaryVertex();
  else
    Printf("GetVertex: ERROR: Invalid second argument %d", mode);
  
  if(!vertex) return 0;
  
  // check Ncontributors
  if(vertex->GetNContributors() <= 0) return 0;
  
  // check resolution
  Double_t zRes = vertex->GetZRes();
  if(zRes == 0) return 0;
  
  //check position
  if(TMath::Abs(vertex->GetXv()) > gVxMax) return 0;
  if(TMath::Abs(vertex->GetYv()) > gVyMax) return 0;
  if(TMath::Abs(vertex->GetZv()) > gVzMax) return 0;
  
  return vertex;
}




