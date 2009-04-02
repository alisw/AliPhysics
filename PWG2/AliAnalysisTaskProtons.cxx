#include "TStyle.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TList.h"
#include "TFile.h"
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
#include "AliVertexerTracks.h"
#include "AliESDVertex.h"

#include "AliProtonAnalysis.h"
#include "AliProtonAnalysisBase.h"
#include "AliAnalysisTaskProtons.h"

// Analysis task to run the \bar{p}/p analysis
// Author: Panos Cristakoglou

ClassImp(AliAnalysisTaskProtons)
  
//________________________________________________________________________ 
AliAnalysisTaskProtons::AliAnalysisTaskProtons()
  : AliAnalysisTask(), fESD(0), fAOD(0), fMC(0),
    fList(0), fProtonAnalysis(0) {
  //Dummy constructor
  
}

//________________________________________________________________________
AliAnalysisTaskProtons::AliAnalysisTaskProtons(const char *name) 
  : AliAnalysisTask(name, ""), fESD(0), fAOD(0), fMC(0),
    fList(0), fProtonAnalysis(0) {
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
  TString gAnalysisLevel = dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetAnalysisLevel(); 

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    if(gAnalysisLevel == "ESD") {
      AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
     
      if (!esdH) {
	Printf("ERROR: Could not get ESDInputHandler");
      } else
	fESD = esdH->GetEvent();
    }
    else if(gAnalysisLevel == "AOD") {
     AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
     
     if (!aodH) {
       Printf("ERROR: Could not get AODInputHandler");
     } else
       fAOD = aodH->GetEvent();
    }
    else if(gAnalysisLevel == "MC") {
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
  // Create output objects
  // Called once
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
  TString gAnalysisLevel = dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetAnalysisLevel(); 
  
  if(gAnalysisLevel == "ESD") {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }

    if(dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->IsEventTriggered(fESD,dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetTriggerMode())) {
      const AliESDVertex *vertex = dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVertex(fESD,dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetAnalysisMode(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVxMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVyMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVzMax());
      if(vertex) {
	Printf("Proton ESD analysis task: There are %d tracks in this event", fESD->GetNumberOfTracks());
	fProtonAnalysis->Analyze(fESD,vertex);
      }//reconstructed vertex
    }//triggered event
  }//ESD analysis              
  
  else if(gAnalysisLevel == "AOD") {
    if (!fAOD) {
      Printf("ERROR: fAOD not available");
      return;
    }
    
    Printf("Proton AOD analysis task: There are %d tracks in this event", fAOD->GetNumberOfTracks());
    fProtonAnalysis->Analyze(fAOD);
  }//AOD analysis
  
  else if(gAnalysisLevel == "MC") {
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
    fProtonAnalysis->Analyze(stack,kFALSE);//kTRUE in case of inclusive measurement
  }//MC analysis                      

  // Post output data.
  PostData(0, fList);
}      

//________________________________________________________________________
void AliAnalysisTaskProtons::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
  gStyle->SetPalette(1,0);

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

  TCanvas *c2 = dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetListOfCuts();
  TFile *flocal = TFile::Open("ListOfCuts.root","recreate");
  c2->Write();
  flocal->Close();
}

