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
//#include "AliESDtrackCuts.h"
#include "AliCentrality.h"

#include "AliEbyEFluctuationAnalysisTaskTrain.h"

// Event by event charge fluctuation analysis
// Authors: Satyajit Jena and Panos Cristakoglou
// 

ClassImp(AliEbyEFluctuationAnalysisTaskTrain)

//________________________________________________________________________
AliEbyEFluctuationAnalysisTaskTrain::AliEbyEFluctuationAnalysisTaskTrain(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fOutputList(0), 
    fHistEventStats(0), fHistCentrality(0),
    fHistNMult(0), fHistNPlusNMinus(0),
    fHistNMultMC(0), fHistNPlusNMinusMC(0),
    fAnalysisType("ESD"), fAnalysisMode("TPC"),
    fCentralityEstimator("V0M"), 
    fCentralityPercentileMin(0.), fCentralityPercentileMax(5.),
    fVxMax(3.), fVyMax(3.), fVzMax(10.),
    fMinNumberOfTPCClusters(80), fMaxChi2PerTPCCluster(4.0), 
    fMaxDCAxy(3.0), fMaxDCAz(3.0) {
  // Constructor
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliEbyEFluctuationAnalysisTaskTrain::UserCreateOutputObjects() {
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
			       20,0.5,20.5);
    fOutputList->Add(fHistCentrality);
    
    TString histName;
    histName = "fHistNMult";
    fHistNMult = new TH1F(histName.Data(), 
			  ";N_{mult.}",
			  500,0,3000);
    fOutputList->Add(fHistNMult);

    histName = "fHistNPlusNMinus";
    fHistNPlusNMinus = new TH2F(histName.Data(), 
				";N_{+};N_{-}",
				2000,0.5,2000.5,2000,0.5,2000.5);  
    fOutputList->Add(fHistNPlusNMinus);
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

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliEbyEFluctuationAnalysisTaskTrain::UserExec(Option_t *) {
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
      //Printf("Event accepted");

      //Centrality stuff
      AliCentrality *centrality = fESD->GetCentrality();
      Int_t nCentrality = 0;
      
      if(centrality->IsEventInCentralityClass(fCentralityPercentileMin,
					      fCentralityPercentileMax,
					      fCentralityEstimator.Data())) {
	if(fAnalysisMode = "TPC") {
	  const AliESDVertex *vertex = fESD->GetPrimaryVertexTPC();
	  if(vertex) {
	    if(vertex->GetNContributors() > 0) {
	      if(vertex->GetZRes() != 0) {
		fHistEventStats->Fill(3); //events with a proper vertex
		if(TMath::Abs(vertex->GetXv()) < fVxMax) {
		  if(TMath::Abs(vertex->GetYv()) < fVyMax) {
		    if(TMath::Abs(vertex->GetZv()) < fVzMax) {
		      fHistEventStats->Fill(4); //analyzed events
		      
		      Printf("Centrality percentile: %lf - Centrality: %d - Total tracks: %d",
			     centrality->GetCentralityPercentile(fCentralityEstimator.Data()),
			     nCentrality,fESD->GetNumberOfTracks());
		      
		      Int_t nAcceptedTracks = 0;
		      //TObjArray *gTrackArray = 0;
		      //if(fESDtrackCuts)
		      //gTrackArray = fESDtrackCuts->GetAcceptedTracks(fESD,kTRUE);
		      //if(gTrackArray) {
		      //nAcceptedTracks = gTrackArray->GetEntries();
		      nAcceptedTracks = fESD->GetNumberOfTracks();
		      AliESDtrack* track = 0;

		      Float_t dcaXY = 0.0, dcaZ = 0.0;

		      // Track loop to fill a pT spectrum
		      for (Int_t iTracks = 0; iTracks < nAcceptedTracks; iTracks++) {
			//track = dynamic_cast<AliESDtrack *>(gTrackArray->At(iTracks));
			track = dynamic_cast<AliESDtrack *>(fESD->GetTrack(iTracks));
			if (!track) {
			  printf("ERROR: Could not receive track %d\n", iTracks);
			  continue;
			}

			AliESDtrack *tpcOnlyTrack = new AliESDtrack();
			
			// only true if we have a tpc track
			if (!track->FillTPCOnlyTrack(*tpcOnlyTrack)) {
			  delete tpcOnlyTrack;
			  return;
			}
			
			tpcOnlyTrack->GetImpactParameters(dcaXY,dcaZ);
			
			//==============================================================//
			Int_t nClustersITS = track->GetITSclusters(0x0);
			Int_t nClustersTPC = track->GetTPCclusters(0x0);

			Float_t chi2PerClusterITS = -1;
			if (nClustersITS!=0)
			  chi2PerClusterITS = track->GetITSchi2()/Float_t(nClustersITS);
			Float_t chi2PerClusterTPC = -1;
			if (nClustersTPC!=0)
			  chi2PerClusterTPC = track->GetTPCchi2()/Float_t(nClustersTPC);
			
			//Printf("TPC clusters: %d - %f === dca: %lf %lf",nClustersTPC, chi2PerClusterTPC,dcaXY,dcaZ);

			if(nClustersTPC < fMinNumberOfTPCClusters) continue;
			if(chi2PerClusterTPC > fMaxChi2PerTPCCluster) continue;
			if(TMath::Abs(dcaXY) > fMaxDCAxy) continue;
			if(TMath::Abs(dcaZ) > fMaxDCAz) continue;
			//==============================================================//

			//if(fESDtrackCuts) {  
			//AliESDtrack *tpcOnlyTrack = fESDtrackCuts->GetTPCOnlyTrack(fESD,iTracks);  
			//if(!tpcOnlyTrack) continue;  
			    
			//if(!fESDtrackCuts->AcceptTrack(tpcOnlyTrack)) continue;  
			//}//ESD track cuts object

			Short_t gCharge = tpcOnlyTrack->Charge();
			if(gCharge > 0) nPlus += 1;
			if(gCharge < 0) nMinus += 1;
		      }//track loop
		      //}//TObjArray valid object
		      //if((nCentrality >= 1)&&(nCentrality <= 20)) {
		      cout<<"=========NPlus: "<<nPlus<<"=========NMinus: "<<nMinus<<endl;
		      //Printf("===NPlus: %d - NMinus: %d===",nPlus,nMinus);
		      fHistCentrality->Fill(nCentrality);
		      fHistNPlusNMinus->Fill(nPlus,nMinus);
		      fHistNMult->Fill(nPlus+nMinus);
		      //}
		    }//Vz cut
		  }//Vy cut
		}//Vx cut
	      }//Vz resolution
	    }//number of contributors
	  }//valid vertex
	}//TPC analysis mode
      }//centrality
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


}      

//________________________________________________________________________
void AliEbyEFluctuationAnalysisTaskTrain::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}
