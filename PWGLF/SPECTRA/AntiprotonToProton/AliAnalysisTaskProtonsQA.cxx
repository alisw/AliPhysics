#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TH1F.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliGenEventHeader.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliESDVertex.h"
#include "AliPhysicsSelection.h"

#include "AliProtonQAAnalysis.h"
#include "AliProtonAnalysisBase.h"
#include "AliAnalysisTaskProtonsQA.h"

//-----------------------------------------------------------------
//                 AliAnalysisTakProtonsQA class
//   This is the task to run the \bar{p}/p QA analysis
//   Origin: Panos Christakoglou | Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

ClassImp(AliAnalysisTaskProtonsQA)
  
//________________________________________________________________________ 
AliAnalysisTaskProtonsQA::AliAnalysisTaskProtonsQA()
  : AliAnalysisTask(), fESD(0), fMC(0), fHistEventStats(0),
    fList0(0), fList1(0), fList2(0), fList3(0), 
    fList4(0), fList5(0), fList6(0), fList7(0), fList8(0),
    fProtonQAAnalysis(0) {
  //Dummy constructor
}

//________________________________________________________________________
AliAnalysisTaskProtonsQA::AliAnalysisTaskProtonsQA(const char *name) 
  : AliAnalysisTask(name, ""), fESD(0), fMC(0), fHistEventStats(0),
    fList0(0), fList1(0), fList2(0), fList3(0), 
    fList4(0), fList5(0), fList6(0), fList7(0), fList8(0),
    fProtonQAAnalysis(0) {
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
  DefineOutput(6, TList::Class());
  DefineOutput(7, TList::Class());
  DefineOutput(8, TList::Class());
  DefineOutput(9, TH1F::Class());
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
  TString gCutName[5] = {"Total","Triggered","Offline trigger",
			 "Vertex","Analyzed"};
  fHistEventStats = new TH1F("fHistEventStats",
			     "Event statistics;;N_{events}",
			     5,0.5,5.5);
  for(Int_t i = 1; i <= 5; i++) 
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());

  fList0 = new TList();
  fList0 = fProtonQAAnalysis->GetGlobalQAList();

  fList1 = new TList();
  fList1 = fProtonQAAnalysis->GetPDGList();

  fList2 = new TList();
  fList2 = fProtonQAAnalysis->GetMCProcessesList();

  fList3 = new TList();
  fList3 = fProtonQAAnalysis->GetAcceptedCutList();

  fList4 = new TList();
  fList4 = fProtonQAAnalysis->GetRejectedCutList();

  fList5 = new TList();
  fList5 = fProtonQAAnalysis->GetAcceptedDCAList();

  fList6 = new TList();
  fList6 = fProtonQAAnalysis->GetEfficiencyQAList();

  fList7 = new TList();
  fList7 = fProtonQAAnalysis->GetVertexQAList();

  fList8 = new TList();
  fList8 = fProtonQAAnalysis->GetCutEfficiencyList();
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
  
  AliGenEventHeader *header = fMC->GenEventHeader();
  if (!header) {
     Printf("ERROR: Could not retrieve the header");
     return;
  }

  AliStack* stack = fMC->Stack();
  if (!stack) {
    Printf("ERROR: Could not retrieve the stack");
    return;
  }

  fHistEventStats->Fill(1);
  if(dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->IsOnlineTriggerUsed()) {
    //online trigger
    if(dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->IsEventTriggered(fESD,dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetTriggerMode())) {
      fHistEventStats->Fill(2);
      
      //offline trigger
      if(dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->IsOfflineTriggerUsed()) {
	AliPhysicsSelection *gPhysicselection = dynamic_cast<AliPhysicsSelection *>(dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetPhysicsSelectionObject());
	if(gPhysicselection->IsCollisionCandidate(fESD)) {
	  fHistEventStats->Fill(3);
	  
	  fProtonQAAnalysis->RunVertexQA(header,
					 fESD);
	  const AliESDVertex *vertex = dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetVertex(fESD,dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetAnalysisMode(),dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetVxMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetVyMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetVzMax());
	  fHistEventStats->Fill(4);
	  
	  if(vertex) {
	    fHistEventStats->Fill(5);
	    fProtonQAAnalysis->RunQAAnalysis(stack, fESD, vertex);
	    fProtonQAAnalysis->RunMCAnalysis(stack);
	    fProtonQAAnalysis->RunPIDEfficiencyAnalysis(stack, fESD, vertex);
	    fProtonQAAnalysis->RunReconstructionEfficiencyAnalysis(fMC,fESD,vertex);
	    fProtonQAAnalysis->RunCutEfficiencyAnalysis(stack, fESD, vertex);
	  }//accepted vertex
	}//offline trigger
      }//offline trigger used
      else {
	fHistEventStats->Fill(3);
	
	fProtonQAAnalysis->RunVertexQA(header,
				       fESD);
	const AliESDVertex *vertex = dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetVertex(fESD,dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetAnalysisMode(),dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetVxMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetVyMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetVzMax());
	fHistEventStats->Fill(4);
	
	if(vertex) {
	  fHistEventStats->Fill(5);
	  fProtonQAAnalysis->RunQAAnalysis(stack, fESD, vertex);
	  fProtonQAAnalysis->RunMCAnalysis(stack);
	  fProtonQAAnalysis->RunPIDEfficiencyAnalysis(stack, fESD, vertex);
	  fProtonQAAnalysis->RunReconstructionEfficiencyAnalysis(fMC,fESD,vertex);
	  fProtonQAAnalysis->RunCutEfficiencyAnalysis(stack, fESD, vertex);
	}//accepted vertex
      }//offline trigger not used
    }//triggered event - online
  }//online trigger used
  else {
    fHistEventStats->Fill(2);
    
    //offline trigger
    if(dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->IsOfflineTriggerUsed()) {
      AliPhysicsSelection *gPhysicselection = dynamic_cast<AliPhysicsSelection *>(dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetPhysicsSelectionObject());
      if(gPhysicselection->IsCollisionCandidate(fESD)) {
	fHistEventStats->Fill(3);
	
	fProtonQAAnalysis->RunVertexQA(header,
				       fESD);
	const AliESDVertex *vertex = dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetVertex(fESD,dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetAnalysisMode(),dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetVxMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetVyMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetVzMax());
	fHistEventStats->Fill(4);
	
	if(vertex) {
	  fHistEventStats->Fill(5);
	  fProtonQAAnalysis->RunQAAnalysis(stack, fESD, vertex);
	  fProtonQAAnalysis->RunMCAnalysis(stack);
	  fProtonQAAnalysis->RunPIDEfficiencyAnalysis(stack, fESD, vertex);
	  fProtonQAAnalysis->RunReconstructionEfficiencyAnalysis(fMC,fESD,vertex);
	  fProtonQAAnalysis->RunCutEfficiencyAnalysis(stack, fESD, vertex);
	}//accepted vertex
      }//offline trigger
    }//offline trigger used
    else {
      fHistEventStats->Fill(3);
      
      fProtonQAAnalysis->RunVertexQA(header,
				     fESD);
      const AliESDVertex *vertex = dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetVertex(fESD,dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetAnalysisMode(),dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetVxMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetVyMax(),dynamic_cast<AliProtonAnalysisBase*>(fProtonQAAnalysis->GetProtonAnalysisBaseObject())->GetVzMax());
      fHistEventStats->Fill(4);
      
      if(vertex) {
	fHistEventStats->Fill(5);
	fProtonQAAnalysis->RunQAAnalysis(stack, fESD, vertex);
	fProtonQAAnalysis->RunMCAnalysis(stack);
	fProtonQAAnalysis->RunPIDEfficiencyAnalysis(stack, fESD, vertex);
	fProtonQAAnalysis->RunReconstructionEfficiencyAnalysis(fMC,fESD,vertex);
	fProtonQAAnalysis->RunCutEfficiencyAnalysis(stack, fESD, vertex);
      }//accepted vertex
    }//offline trigger not used
  }//online trigger not used

  // Post output data.
  PostData(0, fList0);
  PostData(1, fList1);
  PostData(2, fList2);
  PostData(3, fList3);
  PostData(4, fList4);
  PostData(5, fList5);
  PostData(6, fList6);
  PostData(7, fList7);
  PostData(8, fList8);
  PostData(9, fHistEventStats);
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
  fList6 = dynamic_cast<TList*> (GetOutputData(6));
  if (!fList6) {
    Printf("ERROR: fList6 not available");
    return;
  }
  fList7 = dynamic_cast<TList*> (GetOutputData(7));
  if (!fList7) {
    Printf("ERROR: fList7 not available");
    return;
  }
  fList8 = dynamic_cast<TList*> (GetOutputData(8));
  if (!fList8) {
    Printf("ERROR: fList8 not available");
    return;
  }
  fHistEventStats = dynamic_cast<TH1F*> (GetOutputData(9));
  if (!fHistEventStats) {
    Printf("ERROR: fHistEventStats not available");
    return;
  }
}


