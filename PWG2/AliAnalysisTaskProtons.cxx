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
#include "AliVertexerTracks.h"
#include "AliESDVertex.h"

#include "AliProtonAnalysis.h"
#include "AliAnalysisTaskProtons.h"

// Analysis task creating a the 2d y-p_t spectrum of p and antip
// Author: Panos Cristakoglou

ClassImp(AliAnalysisTaskProtons)

//________________________________________________________________________ 
AliAnalysisTaskProtons::AliAnalysisTaskProtons()
  : AliAnalysisTask(), fESD(0), fAOD(0), fMC(0), fAnalysisType("ESD"),
    fList(0), fProtonAnalysis(0),
    fElectronFunction(0), fMuonFunction(0),
    fPionFunction(0), fKaonFunction(0), fProtonFunction(0),
    fFunctionUsed(kFALSE) {
  //Dummy constructor
                                                                                                   
}

//________________________________________________________________________
AliAnalysisTaskProtons::AliAnalysisTaskProtons(const char *name) 
: AliAnalysisTask(name, ""), fESD(0), fAOD(0), fMC(0), fAnalysisType("ESD"),
  fList(0), fProtonAnalysis(0), 
  fElectronFunction(0), fMuonFunction(0),
  fPionFunction(0), fKaonFunction(0), fProtonFunction(0),
  fFunctionUsed(kFALSE),
  fTriggerMode(kMB2), fProtonAnalysisMode(kTPC),
  fVxMax(0), fVyMax(0), fVzMax(0) { 
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
  fProtonAnalysis = new AliProtonAnalysis();

  if(fAnalysisType == "ESD") {
    //Use of TPConly tracks
    if(fProtonAnalysisMode == kTPC) {
      fProtonAnalysis->InitAnalysisHistograms(10, -0.5, 0.5, 16, 0.5, 0.9);
      fProtonAnalysis->UseTPCOnly();
      fProtonAnalysis->SetTPCpid();
      fProtonAnalysis->SetMinTPCClusters(100);
      fProtonAnalysis->SetMaxChi2PerTPCCluster(2.2);
      fProtonAnalysis->SetMaxCov11(0.5);
      fProtonAnalysis->SetMaxCov22(0.5);
      fProtonAnalysis->SetMaxCov33(0.5);
      fProtonAnalysis->SetMaxCov44(0.5);
      fProtonAnalysis->SetMaxCov55(0.5);
      fProtonAnalysis->SetMaxSigmaToVertexTPC(2.0);
      //fProtonAnalysis->SetMaxDCAXYTPC(1.5);
      //fProtonAnalysis->SetMaxDCAZTPC(1.5);
    }
    //Use of HybridTPC tracks
    else if(fProtonAnalysisMode == kHybrid) {
      fProtonAnalysis->InitAnalysisHistograms(10, -0.5, 0.5, 16, 0.5, 0.9);
      fProtonAnalysis->UseHybridTPC();
      fProtonAnalysis->SetTPCpid();
      fProtonAnalysis->SetMinTPCClusters(110);
      fProtonAnalysis->SetMaxChi2PerTPCCluster(2.2);
      fProtonAnalysis->SetMaxCov11(0.5);
      fProtonAnalysis->SetMaxCov22(0.5);
      fProtonAnalysis->SetMaxCov33(0.5);
      fProtonAnalysis->SetMaxCov44(0.5);
      fProtonAnalysis->SetMaxCov55(0.5);
      fProtonAnalysis->SetMaxSigmaToVertex(2.0);
      /*fProtonAnalysis->SetMaxDCAXY(1.5);
	fProtonAnalysis->SetMaxDCAZ(1.5);*/
      fProtonAnalysis->SetPointOnITSLayer6();
      fProtonAnalysis->SetPointOnITSLayer5();
      //fProtonAnalysis->SetPointOnITSLayer4();
      //fProtonAnalysis->SetPointOnITSLayer3();
      fProtonAnalysis->SetPointOnITSLayer2();
      fProtonAnalysis->SetPointOnITSLayer1();
      fProtonAnalysis->SetMinITSClusters(5);
    }
    //Combined tracking
    else if(fProtonAnalysisMode == kGlobal) {
      fProtonAnalysis->InitAnalysisHistograms(10, -0.5, 0.5, 16, 0.5, 0.9);
      fProtonAnalysis->SetMinTPCClusters(110);
      fProtonAnalysis->SetMaxChi2PerTPCCluster(2.2);
      fProtonAnalysis->SetMaxCov11(0.5);
      fProtonAnalysis->SetMaxCov22(0.5);
      fProtonAnalysis->SetMaxCov33(0.5);
      fProtonAnalysis->SetMaxCov44(0.5);
      fProtonAnalysis->SetMaxCov55(0.5);
      fProtonAnalysis->SetMaxSigmaToVertex(2.0);
      //fProtonAnalysis->SetMaxDCAXY(2.0);
      //fProtonAnalysis->SetMaxDCAZ(2.0);
      fProtonAnalysis->SetTPCRefit();
      fProtonAnalysis->SetPointOnITSLayer1();
      fProtonAnalysis->SetPointOnITSLayer2();
      //fProtonAnalysis->SetPointOnITSLayer3();
      //fProtonAnalysis->SetPointOnITSLayer4();
      fProtonAnalysis->SetPointOnITSLayer5();
      fProtonAnalysis->SetPointOnITSLayer6();
      fProtonAnalysis->SetMinITSClusters(5);
      fProtonAnalysis->SetITSRefit();
      fProtonAnalysis->SetESDpid();
    }

    if(fFunctionUsed)
      fProtonAnalysis->SetPriorProbabilityFunctions(fElectronFunction,
					      fMuonFunction,
					      fPionFunction,
					      fKaonFunction,
					      fProtonFunction);
    else
      fProtonAnalysis->SetPriorProbabilities(partFrac);
  }//ESD analysis

  fList = new TList();
  fList->Add(fProtonAnalysis->GetProtonYPtHistogram());
  fList->Add(fProtonAnalysis->GetAntiProtonYPtHistogram());
  fList->Add(fProtonAnalysis->GetEventHistogram());
  fList->Add(fProtonAnalysis->GetProtonContainer());
  fList->Add(fProtonAnalysis->GetAntiProtonContainer());
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

    if(IsEventTriggered(fESD,fTriggerMode)) {
      const AliESDVertex *vertex = GetVertex(fESD,fProtonAnalysisMode,
					     fVxMax,fVyMax,fVzMax);
      if(vertex) {
	Printf("Proton ESD analysis task: There are %d tracks in this event", fESD->GetNumberOfTracks());
	fProtonAnalysis->Analyze(fESD,vertex);
      }//reconstructed vertex
    }//triggered event
  }//ESD analysis              
  
  else if(fAnalysisType == "AOD") {
    if (!fAOD) {
      Printf("ERROR: fAOD not available");
      return;
    }
    
    Printf("Proton AOD analysis task: There are %d tracks in this event", fAOD->GetNumberOfTracks());
    fProtonAnalysis->Analyze(fAOD);
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
    fProtonAnalysis->Analyze(stack);
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

//________________________________________________________________________
Bool_t AliAnalysisTaskProtons::IsEventTriggered(const AliESDEvent *esd, 
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
const AliESDVertex* AliAnalysisTaskProtons::GetVertex(AliESDEvent* esd, 
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


