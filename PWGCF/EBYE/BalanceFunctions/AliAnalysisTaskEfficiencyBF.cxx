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

#include "AliAnalysisTaskEfficiencyBF.h"

// ---------------------------------------------------------------------
//
// Task for calculating the efficiency of the Balance Function 
// for single particles and pairs
//
// Authors: Panos Christakoglou, Michael Weber
// 
// ---------------------------------------------------------------------

ClassImp(AliAnalysisTaskEfficiencyBF)

//________________________________________________________________________
AliAnalysisTaskEfficiencyBF::AliAnalysisTaskEfficiencyBF(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fQAList(0), fOutputList(0), 
  fHistEventStats(0), fHistCentrality(0), fHistNMult(0), 
  fHistGeneratedEtaPtPlus(0), fHistFindableEtaPtPlus(0), 
  fHistReconstructedEtaPtPlus(0), fHistSurvivedEtaPtPlus(0),
  fHistGeneratedEtaPtMinus(0), fHistFindableEtaPtMinus(0), 
  fHistReconstructedEtaPtMinus(0), fHistSurvivedEtaPtMinus(0),
  fHistGeneratedEtaPtPlusControl(0), fHistFindableEtaPtPlusControl(0), 
  fHistReconstructedEtaPtPlusControl(0), fHistSurvivedEtaPtPlusControl(0),
  fHistGeneratedEtaPtMinusControl(0), fHistFindableEtaPtMinusControl(0), 
  fHistReconstructedEtaPtMinusControl(0), fHistSurvivedEtaPtMinusControl(0),
  fHistGeneratedEtaPtPlusPlus(0), fHistFindableEtaPtPlusPlus(0), 
  fHistReconstructedEtaPtPlusPlus(0), fHistSurvivedEtaPtPlusPlus(0),
  fHistGeneratedEtaPtMinusMinus(0), fHistFindableEtaPtMinusMinus(0), 
  fHistReconstructedEtaPtMinusMinus(0), fHistSurvivedEtaPtMinusMinus(0),
  fHistGeneratedEtaPtPlusMinus(0), fHistFindableEtaPtPlusMinus(0), 
  fHistReconstructedEtaPtPlusMinus(0), fHistSurvivedEtaPtPlusMinus(0),
  fESDtrackCuts(0), fAnalysisMode(0), 
  fCentralityEstimator("V0M"), fCentralityPercentileMin(0.0), fCentralityPercentileMax(5.0), 
  fVxMax(3.0), fVyMax(3.0), fVzMax(10.), 
  fMinNumberOfTPCClusters(80), fMaxChi2PerTPCCluster(4.0), fMaxDCAxy(3.0), fMaxDCAz(3.0),
  fMinPt(0.3), fMaxPt(1.5), fMaxEta(0.8){  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskEfficiencyBF::UserCreateOutputObjects() {
  // Create histograms
  // Called once

  fQAList = new TList();
  fQAList->SetName("QAList");
  fQAList->SetOwner();

  fOutputList = new TList();
  fOutputList->SetName("OutputList");
  fOutputList->SetOwner();

  //Event stats.
  TString gCutName[4] = {"Total","Offline trigger",
                         "Vertex","Analyzed"};
  fHistEventStats = new TH1F("fHistEventStats",
                             "Event statistics;;N_{events}",
                             4,0.5,4.5);
  for(Int_t i = 1; i <= 4; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fQAList->Add(fHistEventStats);

  //ESD analysis
  fHistCentrality = new TH1F("fHistCentrality",";Centrality bin;Events",
			     20,0.5,20.5);
  fQAList->Add(fHistCentrality);
  
  //multiplicity (good MC tracks)
  TString histName;
  histName = "fHistNMult";
  fHistNMult = new TH1F(histName.Data(), 
			";N_{mult.}",
			200,0,20000);
  fQAList->Add(fHistNMult);
  
  //eta vs pt for MC positives
  fHistGeneratedEtaPtPlus = new TH2F("fHistGeneratedEtaPtPlus",
				     "Generated positive primaries;#eta;p_{T} (GeV/c)",
				     40,-1.0,1.0,49,0.1,5.0);
  fOutputList->Add(fHistGeneratedEtaPtPlus);
  fHistFindableEtaPtPlus = new TH2F("fHistFindableEtaPtPlus",
				     "Findable positive primaries;#eta;p_{T} (GeV/c)",
				     40,-1.0,1.0,49,0.1,5.0);
  fOutputList->Add(fHistFindableEtaPtPlus);
  fHistReconstructedEtaPtPlus = new TH2F("fHistReconstructedEtaPtPlus",
				     "Reconstructed positive primaries;#eta;p_{T} (GeV/c)",
				     40,-1.0,1.0,49,0.1,5.0);
  fOutputList->Add(fHistReconstructedEtaPtPlus);
  fHistSurvivedEtaPtPlus = new TH2F("fHistSurvivedEtaPtPlus",
				     "Survived positive primaries;#eta;p_{T} (GeV/c)",
				     40,-1.0,1.0,49,0.1,5.0);
  fOutputList->Add(fHistSurvivedEtaPtPlus);

  //eta vs pt for MC negatives
  fHistGeneratedEtaPtMinus = new TH2F("fHistGeneratedEtaPtMinus",
				     "Generated positive primaries;#eta;p_{T} (GeV/c)",
				     40,-1.0,1.0,49,0.1,5.0);
  fOutputList->Add(fHistGeneratedEtaPtMinus);
  fHistFindableEtaPtMinus = new TH2F("fHistFindableEtaPtMinus",
				     "Findable positive primaries;#eta;p_{T} (GeV/c)",
				     40,-1.0,1.0,49,0.1,5.0);
  fOutputList->Add(fHistFindableEtaPtMinus);
  fHistReconstructedEtaPtMinus = new TH2F("fHistReconstructedEtaPtMinus",
				     "Reconstructed positive primaries;#eta;p_{T} (GeV/c)",
				     40,-1.0,1.0,49,0.1,5.0);
  fOutputList->Add(fHistReconstructedEtaPtMinus);
  fHistSurvivedEtaPtMinus = new TH2F("fHistSurvivedEtaPtMinus",
				     "Survived positive primaries;#eta;p_{T} (GeV/c)",
				     40,-1.0,1.0,49,0.1,5.0);
  fOutputList->Add(fHistSurvivedEtaPtMinus);

  //eta vs pt for MC positives (control)
  fHistGeneratedEtaPtPlusControl = new TH2F("fHistGeneratedEtaPtPlusControl",
				     "Generated positive primaries;#eta;p_{T} (GeV/c)",
				     40,-1.0,1.0,49,0.1,5.0);
  fOutputList->Add(fHistGeneratedEtaPtPlusControl);
  fHistFindableEtaPtPlusControl = new TH2F("fHistFindableEtaPtPlusControl",
				     "Findable positive primaries;#eta;p_{T} (GeV/c)",
				     40,-1.0,1.0,49,0.1,5.0);
  fOutputList->Add(fHistFindableEtaPtPlusControl);
  fHistReconstructedEtaPtPlusControl = new TH2F("fHistReconstructedEtaPtPlusControl",
				     "Reconstructed positive primaries;#eta;p_{T} (GeV/c)",
				     40,-1.0,1.0,49,0.1,5.0);
  fOutputList->Add(fHistReconstructedEtaPtPlusControl);
  fHistSurvivedEtaPtPlusControl = new TH2F("fHistSurvivedEtaPtPlusControl",
				     "Survived positive primaries;#eta;p_{T} (GeV/c)",
				     40,-1.0,1.0,49,0.1,5.0);
  fOutputList->Add(fHistSurvivedEtaPtPlusControl);

  //eta vs pt for MC negatives (control)
  fHistGeneratedEtaPtMinusControl = new TH2F("fHistGeneratedEtaPtMinusControl",
				     "Generated positive primaries;#eta;p_{T} (GeV/c)",
				     40,-1.0,1.0,49,0.1,5.0);
  fOutputList->Add(fHistGeneratedEtaPtMinusControl);
  fHistFindableEtaPtMinusControl = new TH2F("fHistFindableEtaPtMinusControl",
				     "Findable positive primaries;#eta;p_{T} (GeV/c)",
				     40,-1.0,1.0,49,0.1,5.0);
  fOutputList->Add(fHistFindableEtaPtMinusControl);
  fHistReconstructedEtaPtMinusControl = new TH2F("fHistReconstructedEtaPtMinusControl",
				     "Reconstructed positive primaries;#eta;p_{T} (GeV/c)",
				     40,-1.0,1.0,49,0.1,5.0);
  fOutputList->Add(fHistReconstructedEtaPtMinusControl);
  fHistSurvivedEtaPtMinusControl = new TH2F("fHistSurvivedEtaPtMinusControl",
				     "Survived positive primaries;#eta;p_{T} (GeV/c)",
				     40,-1.0,1.0,49,0.1,5.0);
  fOutputList->Add(fHistSurvivedEtaPtMinusControl);

  //eta vs pt for MC ++
  fHistGeneratedEtaPtPlusPlus = new TH2F("fHistGeneratedEtaPtPlusPlus",
				     "Generated ++ primaries;#Delta#eta;p_{T} (GeV/c)",
				     40,0.0,2.0,49,0.1,5.0);
  fOutputList->Add(fHistGeneratedEtaPtPlusPlus);
  fHistFindableEtaPtPlusPlus = new TH2F("fHistFindableEtaPtPlusPlus",
				     "Findable ++ primaries;#Delta#eta;p_{T} (GeV/c)",
				     40,0.0,2.0,49,0.1,5.0);
  fOutputList->Add(fHistFindableEtaPtPlusPlus);
  fHistReconstructedEtaPtPlusPlus = new TH2F("fHistReconstructedEtaPtPlusPlus",
				     "Reconstructed ++ primaries;#Delta#eta;p_{T} (GeV/c)",
				     40,0.0,2.0,49,0.1,5.0);
  fOutputList->Add(fHistReconstructedEtaPtPlusPlus);
  fHistSurvivedEtaPtPlusPlus = new TH2F("fHistSurvivedEtaPtPlusPlus",
				     "Survived ++ primaries;#Delta#eta;p_{T} (GeV/c)",
				     40,0.0,2.0,49,0.1,5.0);
  fOutputList->Add(fHistSurvivedEtaPtPlusPlus);

  //eta vs pt for MC --
  fHistGeneratedEtaPtMinusMinus = new TH2F("fHistGeneratedEtaPtMinusMinus",
				     "Generated -- primaries;#Delta#eta;p_{T} (GeV/c)",
				     40,0.0,2.0,49,0.1,5.0);
  fOutputList->Add(fHistGeneratedEtaPtMinusMinus);
  fHistFindableEtaPtMinusMinus = new TH2F("fHistFindableEtaPtMinusMinus",
				     "Findable -- primaries;#Delta#eta;p_{T} (GeV/c)",
				     40,0.0,2.0,49,0.1,5.0);
  fOutputList->Add(fHistFindableEtaPtMinusMinus);
  fHistReconstructedEtaPtMinusMinus = new TH2F("fHistReconstructedEtaPtMinusMinus",
				     "Reconstructed -- primaries;#Delta#eta;p_{T} (GeV/c)",
				     40,0.0,2.0,49,0.1,5.0);
  fOutputList->Add(fHistReconstructedEtaPtMinusMinus);
  fHistSurvivedEtaPtMinusMinus = new TH2F("fHistSurvivedEtaPtMinusMinus",
				     "Survived -- primaries;#Delta#eta;p_{T} (GeV/c)",
				     40,0.0,2.0,49,0.1,5.0);
  fOutputList->Add(fHistSurvivedEtaPtMinusMinus);

  //eta vs pt for MC +-
  fHistGeneratedEtaPtPlusMinus = new TH2F("fHistGeneratedEtaPtPlusMinus",
				     "Generated +- primaries;#Delta#eta;p_{T} (GeV/c)",
				     40,0.0,2.0,49,0.1,5.0);
  fOutputList->Add(fHistGeneratedEtaPtPlusMinus);
  fHistFindableEtaPtPlusMinus = new TH2F("fHistFindableEtaPtPlusMinus",
				     "Findable +- primaries;#Delta#eta;p_{T} (GeV/c)",
				     40,0.0,2.0,49,0.1,5.0);
  fOutputList->Add(fHistFindableEtaPtPlusMinus);
  fHistReconstructedEtaPtPlusMinus = new TH2F("fHistReconstructedEtaPtPlusMinus",
				     "Reconstructed +- primaries;#Delta#eta;p_{T} (GeV/c)",
				     40,0.0,2.0,49,0.1,5.0);
  fOutputList->Add(fHistReconstructedEtaPtPlusMinus);
  fHistSurvivedEtaPtPlusMinus = new TH2F("fHistSurvivedEtaPtPlusMinus",
				     "Survived +- primaries;#Delta#eta;p_{T} (GeV/c)",
				     40,0.0,2.0,49,0.1,5.0);
  fOutputList->Add(fHistSurvivedEtaPtPlusMinus);

  fQAList->Print();
  fOutputList->Print();

  PostData(1, fQAList);
  PostData(2, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskEfficiencyBF::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  // Post output data.
  //ESD analysis
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }
  
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

  // arrays for 2 particle histograms
  Int_t nMCLabelCounter         = 0;
  const Int_t maxMCLabelCounter = 20000;

  Double_t pt[maxMCLabelCounter];
  Double_t eta[maxMCLabelCounter];
  Int_t level[maxMCLabelCounter];
  Int_t charge[maxMCLabelCounter];


  //AliInfo(Form("%d %d",mcEvent->GetNumberOfTracks(),fESD->GetNumberOfTracks()));
  fHistEventStats->Fill(1); //all events
    
  //Centrality stuff
  AliCentrality *centrality = fESD->GetCentrality();
  Int_t nCentrality = 0;
  nCentrality = (Int_t)(centrality->GetCentralityPercentile(fCentralityEstimator.Data())/10.);

  //Printf("Centrality: %lf",centrality->GetCentralityPercentile(fCentralityEstimator.Data()));

  if(centrality->IsEventInCentralityClass(fCentralityPercentileMin,
					  fCentralityPercentileMax,
					  fCentralityEstimator.Data())) {
    fHistEventStats->Fill(2); //triggered + centrality
    fHistCentrality->Fill(nCentrality+1);

    //Printf("Centrality selection: %lf - %lf",fCentralityPercentileMin,fCentralityPercentileMax);
  
    if(fAnalysisMode.CompareTo("TPC") == 0 ) {
      const AliESDVertex *vertex = fESD->GetPrimaryVertexTPC();
      if(vertex) {
	if(vertex->GetNContributors() > 0) {
	  if(vertex->GetZRes() != 0) {
	    fHistEventStats->Fill(3); //events with a proper vertex
	    if(TMath::Abs(vertex->GetXv()) < fVxMax) {
	      if(TMath::Abs(vertex->GetYv()) < fVyMax) {
		if(TMath::Abs(vertex->GetZv()) < fVzMax) {
		  fHistEventStats->Fill(4); //analyzed events
		  
		  Int_t nMCParticles = mcEvent->GetNumberOfTracks();
		  TArrayI labelMCArray(nMCParticles);
		  
		  for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) {
		    AliMCParticle *mcTrack = (AliMCParticle*) mcEvent->GetTrack(iTracks);
		    if (!mcTrack) {
		      Printf("ERROR: Could not receive track %d (mc loop)", iTracks);
		      continue;
		    }
		    
		    //exclude particles generated out of the acceptance
		    Double_t vz = mcTrack->Zv();
		    if (TMath::Abs(vz) > 50.) continue;
		    
		    //acceptance
		    if(TMath::Abs(mcTrack->Eta()) > 2.5) 
		      continue;
		    if((mcTrack->Pt() > 5.0)||(mcTrack->Pt() < 0.1)) 
		      continue;
		    
		    TParticle* particle = mcTrack->Particle();
		    if(!particle) continue;
		    if(!stack->IsPhysicalPrimary(iTracks)) continue;

		    if(iTracks <= stack->GetNprimary()) {		      
		      Short_t gMCCharge = mcTrack->Charge();
		      
		      if(gMCCharge > 0)
			fHistGeneratedEtaPtPlus->Fill(particle->Eta(),
						      particle->Pt());
		      else if(gMCCharge < 0)
			fHistGeneratedEtaPtMinus->Fill(particle->Eta(),
						       particle->Pt());

		      
		      // findable tracks --> DOES NOT WORK????
		      // Loop over Track References
		      Bool_t labelTPC = kTRUE;
		      AliTrackReference* trackRef = 0;

		      for (Int_t iTrackRef = 0; iTrackRef < mcTrack->GetNumberOfTrackReferences(); iTrackRef++) {
			trackRef = mcTrack->GetTrackReference(iTrackRef);
			if(trackRef) {
			  Int_t detectorId = trackRef->DetectorId();
			  if (detectorId == AliTrackReference::kTPC) {
			    labelTPC = kTRUE;
			    break;
			  }
			}
		      }//loop over track references

		      if(labelTPC) {
			labelMCArray.AddAt(iTracks,nMCLabelCounter);

			if(nMCLabelCounter >= maxMCLabelCounter){
			  AliWarning(Form("MC Label Counter > Limit (%d) --> stop loop here",maxMCLabelCounter));
			  break;
			}

			// fill the arrays for 2 particle analysis
			eta[nMCLabelCounter]    = particle->Eta();
			pt[nMCLabelCounter]     = particle->Pt();
			charge[nMCLabelCounter] = gMCCharge;
			// findable = generated in this case!
			level[nMCLabelCounter]  = 1;
			
			nMCLabelCounter += 1;
			
			if(gMCCharge > 0)
			  fHistFindableEtaPtPlus->Fill(particle->Eta(),
						       particle->Pt());
			else if(gMCCharge < 0)
			  fHistFindableEtaPtMinus->Fill(particle->Eta(),
							particle->Pt());
		      }
		    }//primaries
		  }//loop over MC particles
		
		  fHistNMult->Fill(nMCLabelCounter);
		  
		  // not used so far
		  //Float_t dcaXY = 0.0, dcaZ = 0.0;

		  //ESD track loop
		  Int_t nGoodTracks = fESD->GetNumberOfTracks();
		  
		  TArrayI labelArray(nGoodTracks);
		  Int_t labelCounter = 0;
		  for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
		    AliESDtrack* track = fESD->GetTrack(iTracks);
		    //AliESDtrack* track = fESDtrackCuts->GetTPCOnlyTrack(fESD,iTracks);
		    if(!track) continue;

		    AliESDtrack *tpcOnlyTrack = new AliESDtrack();
		    
		    if (!track->FillTPCOnlyTrack(*tpcOnlyTrack)) {
		      delete tpcOnlyTrack;
		      continue;
		    }

		    Int_t label = TMath::Abs(track->GetTPCLabel());
		    if(IsLabelUsed(labelArray,label)) continue;
		    labelArray.AddAt(label,labelCounter);
		    labelCounter += 1;

		    Bool_t iFound = kFALSE;
		    Int_t mcGoods = nMCLabelCounter;
		    for (Int_t k = 0; k < mcGoods; k++) {
		      Int_t mcLabel = labelMCArray.At(k);
		      iFound = kFALSE;
		    			      
		      if (mcLabel != TMath::Abs(label)) continue;
		      if(mcLabel != label) continue;
		      if(label > stack->GetNtrack()) continue;
		      
		      TParticle *particle = stack->Particle(label);
		      if(!particle) continue;
		      
		      //acceptance
		      if(TMath::Abs(particle->Eta()) > 2.5) 
			continue;
		      if((particle->Pt() > 5.0)||(particle->Pt() < 0.1)) 
			continue;
		      if(!stack->IsPhysicalPrimary(label)) continue;
		      
		      if(label <= stack->GetNprimary()) {
			
			// reconstructed
			level[k]  = 2;
			
			Short_t gCharge = track->Charge();
			if(gCharge > 0)
			  fHistReconstructedEtaPtPlus->Fill(particle->Eta(),
							    particle->Pt());
			else if(gCharge < 0)
			  fHistReconstructedEtaPtMinus->Fill(particle->Eta(),
							     particle->Pt());
			
			// track cuts + analysis kinematic cuts
			if(fESDtrackCuts->AcceptTrack(track) && TMath::Abs(track->Eta()) < fMaxEta && track->Pt() > fMinPt && track->Pt() < fMaxPt ){

			  // survived
			  level[k]  = 3;

			  if(gCharge > 0)
			    fHistSurvivedEtaPtPlus->Fill(particle->Eta(),
							 particle->Pt());
			  else if(gCharge < 0)
			    fHistSurvivedEtaPtMinus->Fill(particle->Eta(),
							  particle->Pt());
			}//track cuts
		      }//primary particles
		    }//findable track loop
		  }//ESD track loop
		
		  labelMCArray.Reset();
		  labelArray.Reset();

		}//Vz cut
	      }//Vy cut
	    }//Vx cut
	  }//Vz resolution
	}//number of contributors
      }//valid vertex
    }//TPC analysis mode
  }//centrality  
      


  // Here comes the 2 particle analysis
  // loop over all good MC particles
  for (Int_t i = 0; i < nMCLabelCounter ; i++) {

    // control 1D histograms (charge might be different?)
    if(charge[i] > 0){
      if(level[i] > 0) fHistGeneratedEtaPtPlusControl->Fill(eta[i],pt[i]);
      if(level[i] > 1) fHistReconstructedEtaPtPlusControl->Fill(eta[i],pt[i]);
      if(level[i] > 2) fHistSurvivedEtaPtPlusControl->Fill(eta[i],pt[i]);
    }
    else if(charge[i] < 0){
      if(level[i] > 0) fHistGeneratedEtaPtMinusControl->Fill(eta[i],pt[i]);
      if(level[i] > 1) fHistReconstructedEtaPtMinusControl->Fill(eta[i],pt[i]);
      if(level[i] > 2) fHistSurvivedEtaPtMinusControl->Fill(eta[i],pt[i]);
    }


    for (Int_t j = i+1; j < nMCLabelCounter ; j++) {
      
      if(charge[i] > 0 && charge[j] > 0 ){
	if(level[i] > 0 && level[j] > 0) fHistGeneratedEtaPtPlusPlus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	if(level[i] > 1 && level[j] > 1) fHistReconstructedEtaPtPlusPlus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	if(level[i] > 2 && level[j] > 2) fHistSurvivedEtaPtPlusPlus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
      }
      
      else if(charge[i] < 0 && charge[j] < 0 ){
	if(level[i] > 0 && level[j] > 0) fHistGeneratedEtaPtMinusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	if(level[i] > 1 && level[j] > 1) fHistReconstructedEtaPtMinusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	if(level[i] > 2 && level[j] > 2) fHistSurvivedEtaPtMinusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
      }
      else if((charge[i] > 0 && charge[j] < 0)||(charge[i] < 0 && charge[j] > 0)){
	if(level[i] > 0 && level[j] > 0) fHistGeneratedEtaPtPlusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	if(level[i] > 1 && level[j] > 1) fHistReconstructedEtaPtPlusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
	if(level[i] > 2 && level[j] > 2) fHistSurvivedEtaPtPlusMinus->Fill(TMath::Abs(eta[i]-eta[j]),pt[i]);
      }
    }
  }	  
}
//________________________________________________________________________
void AliAnalysisTaskEfficiencyBF::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}

//____________________________________________________________________//
Bool_t AliAnalysisTaskEfficiencyBF::IsLabelUsed(TArrayI labelArray, Int_t label) {
  //Checks if the label is used already
  Bool_t status = kFALSE;
  for(Int_t i = 0; i < labelArray.GetSize(); i++) {
    if(labelArray.At(i) == label)
      status = kTRUE;
  }

  return status;
}
