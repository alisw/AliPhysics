#include "Riostream.h"
#include "TStyle.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TList.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliESDVertex.h"

#include "AliProtonAnalysis.h"
#include "AliProtonAnalysisBase.h"
#include "AliAnalysisTaskProtons.h"

//-----------------------------------------------------------------
//                 AliAnalysisTakProtons class
//   This is the task to run the \bar{p}/p analysis
//   Origin: Panos Christakoglou | Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

ClassImp(AliAnalysisTaskProtons)
  
//________________________________________________________________________ 
AliAnalysisTaskProtons::AliAnalysisTaskProtons()
  : AliAnalysisTask(), fESD(0), fAOD(0), fMC(0),
    fListAnalysis(0), fListQA(0), fHistEventStats(0), 
  fProtonAnalysis(0) {//, fCutCanvas(0) {
  //Dummy constructor
  
}

//________________________________________________________________________
AliAnalysisTaskProtons::AliAnalysisTaskProtons(const char *name) 
  : AliAnalysisTask(name, ""), fESD(0), fAOD(0), fMC(0),
    fListAnalysis(0), fListQA(0), fHistEventStats(0), 
  fProtonAnalysis(0) {//, fCutCanvas(0) {
  // Constructor
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TCanvas::Class());
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
  char *gCutName[5] = {"Total","Triggered","Offline trigger",
		       "Vertex","Analyzed"};
  fHistEventStats = new TH1F("fHistEventStats",
			     "Event statistics;;N_{events}",
			     5,0.5,5.5);
  for(Int_t i = 1; i <= 5; i++) 
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1]);

  fListAnalysis = new TList();
  fListAnalysis->Add(fProtonAnalysis->GetProtonYPtHistogram());
  fListAnalysis->Add(fProtonAnalysis->GetAntiProtonYPtHistogram());
  fListAnalysis->Add(fProtonAnalysis->GetEventHistogram());
  fListAnalysis->Add(fProtonAnalysis->GetProtonContainer());
  fListAnalysis->Add(fProtonAnalysis->GetAntiProtonContainer());
  fListAnalysis->Add(fHistEventStats);

  //fListQA = new TList();
  //fListQA->SetName("fListQA");
  //fListQA->Add(fProtonAnalysis->GetQAList());
  //fListQA->Add(dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVertexQAList());

  //fCutCanvas = dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetListOfCuts();
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
    
    fHistEventStats->Fill(1);
    if(dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->IsOnlineTriggerUsed()) {
      //online trigger
      if(dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->IsEventTriggered(fESD,dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetTriggerMode())) {
	fHistEventStats->Fill(2);
	AliDebug(1,Form("Fired trigger class: %s",fESD->GetFiredTriggerClasses().Data()));
	
	//offline trigger
	if(dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->IsOfflineTriggerUsed()) {
	  AliPhysicsSelection *gPhysicselection = dynamic_cast<AliPhysicsSelection *>(dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetPhysicsSelectionObject());
	  if(gPhysicselection->IsCollisionCandidate(fESD)) {
	    fHistEventStats->Fill(3);
	    AliDebug(1,Form("Fired trigger class: %s",fESD->GetFiredTriggerClasses().Data()));
	    //Reconstructed vertex
	    const AliESDVertex *vertex = dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVertex(fESD,dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetAnalysisMode(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVxMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVyMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVzMax());
	    fHistEventStats->Fill(4);
	    if(vertex) {
	      AliDebug(1,Form("Proton ESD analysis task: There are %d tracks in this event",fESD->GetNumberOfTracks()));
	      fProtonAnalysis->Analyze(fESD,vertex);
	      fHistEventStats->Fill(5);
	    }//reconstructed vertex
	  }//offline trigger
	}//usage of the offline trigger
	else {
	  fHistEventStats->Fill(3);
	  AliDebug(1,Form("Fired trigger class: %s",fESD->GetFiredTriggerClasses().Data()));
	  //Reconstructed vertex
	  const AliESDVertex *vertex = dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVertex(fESD,dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetAnalysisMode(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVxMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVyMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVzMax());
	  fHistEventStats->Fill(4);
	  if(vertex) {
	    AliDebug(1,Form("Proton ESD analysis task: There are %d tracks in this event",fESD->GetNumberOfTracks()));
	    fProtonAnalysis->Analyze(fESD,vertex);
	    fHistEventStats->Fill(5);
	    }//reconstructed vertex
	}//offline trigger not used
      }//triggered event - online
    }//online trigger used
    else {
      fHistEventStats->Fill(2);
      AliDebug(1,Form("Fired trigger class: %s",fESD->GetFiredTriggerClasses().Data()));
      
      //offline trigger
      if(dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->IsOfflineTriggerUsed()) {
	AliPhysicsSelection *gPhysicselection = dynamic_cast<AliPhysicsSelection *>(dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetPhysicsSelectionObject());
	if(gPhysicselection->IsCollisionCandidate(fESD)) {
	  fHistEventStats->Fill(3);
	  AliDebug(1,Form("Fired trigger class: %s",fESD->GetFiredTriggerClasses().Data()));
	  //Reconstructed vertex
	  const AliESDVertex *vertex = dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVertex(fESD,dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetAnalysisMode(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVxMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVyMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVzMax());
	  fHistEventStats->Fill(4);
	  if(vertex) {
	    AliDebug(1,Form("Proton ESD analysis task: There are %d tracks in this event",fESD->GetNumberOfTracks()));
	    fProtonAnalysis->Analyze(fESD,vertex);
	    fHistEventStats->Fill(5);
	  }//reconstructed vertex
	}//offline trigger
      }//usage of the offline trigger
      else {
	fHistEventStats->Fill(3);
	AliDebug(1,Form("Fired trigger class: %s",fESD->GetFiredTriggerClasses().Data()));
	//Reconstructed vertex
	const AliESDVertex *vertex = dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVertex(fESD,dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetAnalysisMode(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVxMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVyMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetVzMax());
	fHistEventStats->Fill(4);
	if(vertex) {
	  AliDebug(1,Form("Proton ESD analysis task: There are %d tracks in this event",fESD->GetNumberOfTracks()));
	  fProtonAnalysis->Analyze(fESD,vertex);
	  fHistEventStats->Fill(5);
	}//reconstructed vertex
      }//offline trigger not used
    }//online trigger not used
  }//ESD analysis              
  
  else if(gAnalysisLevel == "AOD") {
    if (!fAOD) {
      Printf("ERROR: fAOD not available");
      return;
    }
    AliDebug(1,Form("Proton AOD analysis task: There are %d tracks in this event", fAOD->GetNumberOfTracks()));
    //Printf("Proton AOD analysis task: There are %d tracks in this event", fAOD->GetNumberOfTracks());
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
    AliDebug(1,Form("Proton MC analysis task: There are %d primaries in this event", stack->GetNprimary()));
    //Printf("Proton MC analysis task: There are %d primaries in this event", stack->GetNprimary());
    fProtonAnalysis->Analyze(stack,kFALSE);//kTRUE in case of inclusive measurement
  }//MC analysis                      

  // Post output data.
  PostData(0, fListAnalysis);
  //PostData(1, fListQA);
  //PostData(2, fCutCanvas);
}      

//________________________________________________________________________
void AliAnalysisTaskProtons::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
  gStyle->SetPalette(1,0);

  fListAnalysis = dynamic_cast<TList*> (GetOutputData(0));
  if (!fListAnalysis) {
    Printf("ERROR: fListAnalysis not available");
    return;
  }
   
  TH2F *fHistYPtProtons = (TH2F *)fListAnalysis->At(0);
  TH2F *fHistYPtAntiProtons = (TH2F *)fListAnalysis->At(1);
    
  TCanvas *c1 = new TCanvas("c1","p-\bar{p}",200,0,800,400);
  c1->SetFillColor(10); c1->SetHighLightColor(10); c1->Divide(2,1);

  c1->cd(1)->SetLeftMargin(0.15); c1->cd(1)->SetBottomMargin(0.15);  
  if (fHistYPtProtons) fHistYPtProtons->DrawCopy("colz");
  c1->cd(2)->SetLeftMargin(0.15); c1->cd(2)->SetBottomMargin(0.15);  
  if (fHistYPtAntiProtons) fHistYPtAntiProtons->DrawCopy("colz");

  /*TCanvas *c2 = dynamic_cast<AliProtonAnalysisBase*>(fProtonAnalysis->GetProtonAnalysisBaseObject())->GetListOfCuts();
  TFile *flocal = TFile::Open("ListOfCuts.root","recreate");
  c2->Write();
  flocal->Close();*/
}

