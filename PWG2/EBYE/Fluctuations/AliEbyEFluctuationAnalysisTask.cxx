#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TObjArray.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"

#include "AliEbyEFluctuationAnalysisTask.h"

// Event by event charge fluctuation analysis
// Authors: Satyajit Jena and Panos Cristakoglou
// 

ClassImp(AliEbyEFluctuationAnalysisTask)

//________________________________________________________________________
AliEbyEFluctuationAnalysisTask::AliEbyEFluctuationAnalysisTask(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fOutputList(0), 
    fHistEventStats(0), fHistCentrality(0),
    fHistNMultMC(0), fHistNPlusNMinusMC(0), 
    fESDtrackCuts(0),
    fAnalysisType("ESD"), fAnalysisMode("TPC"),
    fCentralityEstimator("V0M"), fCentralityBins20(kFALSE),
    fVxMax(3.), fVyMax(3.), fVzMax(10.) {
  // Constructor
  
  for(Int_t iBin = 1; iBin <= nCentralityBins; iBin++) {
    fHistNMult[iBin-1] = NULL;
    fHistNPlusNMinus[iBin-1] = NULL;
  }

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliEbyEFluctuationAnalysisTask::UserCreateOutputObjects() {
  // Create histograms
  // Called once

  fOutputList = new TList();
  fOutputList->SetOwner();

  //Event stats.
  TString gCutName[4] = {"Total","Offline trigger",
                         "Vertex","Analyzed"};
  fHistEventStats = new TH1F("fHistEventStats",
                             "Event statistics;;N_{events}",
                             4,0.5,4.5);
  for(Int_t i = 1; i <= 4; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fOutputList->Add(fHistEventStats);

  //ESD analysis
  if(fAnalysisType == "ESD") {
    fHistCentrality = new TH1F("fHistCentrality",";Centrality bin;Events",
			       nCentralityBins,0.5,nCentralityBins+0.5);
    fOutputList->Add(fHistCentrality);
    
    TString histName;
    for(Int_t iBin = 1; iBin <= nCentralityBins; iBin++) {
      histName = "fHistNMult"; histName += "Centrality"; histName += iBin; 
      fHistNMult[iBin-1] = new TH1F(histName.Data(), 
				    ";N_{mult.}",
				    500,0,3000);
      fOutputList->Add(fHistNMult[iBin-1]);
    }
    for(Int_t iBin = 1; iBin <= nCentralityBins; iBin++) {
      histName = "fHistNPlusNMinus"; histName += "Centrality"; histName += iBin; 
      fHistNPlusNMinus[iBin-1] = new TH2F(histName.Data(), 
					  ";N_{+};N_{-}",
					  2000,0.5,2000.5,2000,0.5,2000.5);  
      fOutputList->Add(fHistNPlusNMinus[iBin-1]);
    }
  }//ESD analysis
  else if(fAnalysisType == "MC") {
    TString histName = "fHistNMultMC";
    fHistNMultMC = new TH1F(histName.Data(), 
			    ";N_{mult.}",
			    600,0,6000);
    fOutputList->Add(fHistNMultMC);
    
    histName = "fHistNPlusNMinusMC";
    fHistNPlusNMinusMC = new TH2F(histName.Data(), 
				  ";N_{+};N_{-}",
				  3000,0.5,3000.5,3000,0.5,3000.5);  
    fOutputList->Add(fHistNPlusNMinusMC);
  }//MC analysis
}

//________________________________________________________________________
void AliEbyEFluctuationAnalysisTask::UserExec(Option_t *) {
  // Main loop
  // Called for each event
  
  Int_t nPlus = 0, nMinus = 0;

  // Post output data.
  //ESD analysis
  if(fAnalysisType == "ESD") {
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD) {
      printf("ERROR: fESD not available\n");
      return;
    }
  
    fHistEventStats->Fill(1); //all events
    
    //Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
    //if(isSelected) {
      fHistEventStats->Fill(2); //triggered + centrality
      Printf("Event accepted");

      //Centrality stuff
      AliCentrality *centrality = fESD->GetCentrality();
      Int_t nCentrality = 0;
      if(fCentralityBins20) 
	nCentrality = centrality->GetCentralityClass5(fCentralityEstimator.Data());
      else
	nCentrality = centrality->GetCentralityClass10(fCentralityEstimator.Data());

      if((nCentrality < 0)||(nCentrality > 19)) return;
      
      if(fAnalysisMode == "TPC") {
	const AliESDVertex *vertex = fESD->GetPrimaryVertexTPC();
	if(vertex) {
	  if(vertex->GetNContributors() > 0) {
	    if(vertex->GetZRes() != 0) {
	      fHistEventStats->Fill(3); //events with a proper vertex
	      if(TMath::Abs(vertex->GetXv()) < fVxMax) {
		if(TMath::Abs(vertex->GetYv()) < fVyMax) {
		  if(TMath::Abs(vertex->GetZv()) < fVzMax) {
		    fHistEventStats->Fill(4); //analyzed events
		    
		    //Printf("Centrality percentile: %lf - Centrality: %d - Total tracks: %d",
		    //centrality->GetCentralityPercentile(fCentralityEstimator.Data()),
		    //nCentrality,fESD->GetNumberOfTracks());
		    		    
		    Int_t nAcceptedTracks = 0;
		    TObjArray *gTrackArray = 0;
		    if(fESDtrackCuts)
                      gTrackArray = fESDtrackCuts->GetAcceptedTracks(fESD,kTRUE);
                    if(gTrackArray) {
                      nAcceptedTracks = gTrackArray->GetEntries();

		      AliESDtrack* track = 0;
		      // Track loop to fill a pT spectrum
		      for (Int_t iTracks = 0; iTracks < nAcceptedTracks; iTracks++) {
			track = dynamic_cast<AliESDtrack *>(gTrackArray->At(iTracks));
			if (!track) {
			  printf("ERROR: Could not receive track %d\n", iTracks);
			  continue;
			}
			
			Short_t gCharge = track->Charge();
			if(gCharge > 0) nPlus += 1;
			if(gCharge < 0) nMinus += 1;
		      }//track loop
		    }//TObjArray valid object
		    //if((nCentrality >= 1)&&(nCentrality <= 20)) {
		    
		    fHistCentrality->Fill(nCentrality);
		    fHistNPlusNMinus[nCentrality-1]->Fill(nPlus,nMinus);
		    fHistNMult[nCentrality-1]->Fill(nPlus+nMinus);
		    //}
		  }//Vz cut
		}//Vy cut
	      }//Vx cut
	    }//Vz resolution
	  }//number of contributors
	}//valid vertex
      }//TPC analysis mode
      //}//physics selection
  }//ESD analysis level
  
  //MC analysis
  if(fAnalysisType == "MC") {
    AliMCEvent* mcEvent = MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }
    AliStack* stack = mcEvent->Stack();
    if (!stack) {
      Printf("ERROR: Could not retrieve MC stack");
      return;
    }
    
    fHistEventStats->Fill(1); 
    fHistEventStats->Fill(2); 
    fHistEventStats->Fill(3); 
    fHistEventStats->Fill(4); //analyzed events
    for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) {
      AliVParticle* track = mcEvent->GetTrack(iTracks);
      if (!track) {
	Printf("ERROR: Could not receive track %d (mc loop)", iTracks);
	continue;
      }
      
      if(!stack->IsPhysicalPrimary(iTracks)) continue;
      if((track->Pt() < 0.3) && (track->Pt() > 1.5)) continue;
      if(TMath::Abs(track->Eta()) > 0.8) continue;
      
      Short_t gCharge = track->Charge();
      if(gCharge > 0) nPlus += 1;
      if(gCharge < 0) nMinus += 1;
    }//particle loop
    fHistNPlusNMinusMC->Fill(nPlus,nMinus);
    fHistNMultMC->Fill(nPlus+nMinus);
    
  }//MC analysis level


  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliEbyEFluctuationAnalysisTask::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}
