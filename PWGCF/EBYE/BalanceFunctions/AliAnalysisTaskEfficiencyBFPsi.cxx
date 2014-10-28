#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TObjArray.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliTHn.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliGenHijingEventHeader.h"

#include "AliAnalysisTaskEfficiencyBFPsi.h"

// ---------------------------------------------------------------------
//
// Task for calculating the efficiency of the Balance Function 
// for single particles and pairs
//
// Authors: Panos Christakoglou, Michael Weber
// 
// ---------------------------------------------------------------------

ClassImp(AliAnalysisTaskEfficiencyBFPsi)
  
//________________________________________________________________________
  AliAnalysisTaskEfficiencyBFPsi::AliAnalysisTaskEfficiencyBFPsi(const char *name): 
    AliAnalysisTaskSE(name), fESD(0), fQAList(0), fOutputList(0), 
    fHistEventStats(0), fHistCentrality(0), fHistNMult(0), 
    fHistGeneratedPlus(0), fHistSurvivedPlus(0),
    fHistGeneratedMinus(0), fHistSurvivedMinus(0),
    fHistGeneratedPlusPlus(0), fHistSurvivedPlusPlus(0),
    fHistGeneratedMinusMinus(0), fHistSurvivedMinusMinus(0),
    fHistGeneratedPlusMinus(0), fHistSurvivedPlusMinus(0),
    fHistGeneratedMinusPlus(0), fHistSurvivedMinusPlus(0),
    fESDtrackCuts(0), fAnalysisMode(0), 
    fCentralityEstimator("V0M"), 
    fCentralityPercentileMin(0.0), fCentralityPercentileMax(5.0), 
    fVxMax(3.0), fVyMax(3.0), fVzMax(10.), 
    fMinNumberOfTPCClusters(80), fMaxChi2PerTPCCluster(4.0), 
    fMaxDCAxy(3.0), fMaxDCAz(3.0),
    fMinPt(0.3), fMaxPt(1.5), fMaxEta(0.8), fEtaRangeMax(0.8), 
    fPtRangeMin(0.1), fPtRangeMax(5.0), fPhiRangeMin(0.0),fPhiRangeMax(6.28) {  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskEfficiencyBFPsi::UserCreateOutputObjects() {
  // Create histograms
  // Called once

  // global switch disabling the reference 
  // (to avoid "Replacing existing TH1" if several wagons are created in train)
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

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
  
  //Setting up the AliTHn
  Int_t anaSteps   = 1;       // analysis steps
  Int_t iBinSingle[kVariablesSingle];        // binning for track variables
  Double_t* dBinsSingle[kVariablesSingle];   // bins for track variables  
  TString axisTitleSingle[kVariablesSingle]; // axis titles for track variables
  
  // two particle histograms
  Int_t iBinPair[kVariablesPair];         // binning for track variables
  Double_t* dBinsPair[kVariablesPair];    // bins for track variables  
  TString axisTitlePair[kVariablesPair];  // axis titles for track variables

  //Psi_2: -0.5->0.5 (in plane), 0.5->1.5 (intermediate), 1.5->2.5 (out of plane), 2.5->3.5 (all)
  const Int_t kNPsi2Bins = 4;
  Double_t psi2Bins[kNPsi2Bins+1] = {-0.5,0.5,1.5,2.5,3.5};
  iBinSingle[0]       = kNPsi2Bins;
  dBinsSingle[0]      = psi2Bins;
  axisTitleSingle[0]  = "#phi - #Psi_{2} (a.u.)";
  iBinPair[0]       = kNPsi2Bins;
  dBinsPair[0]      = psi2Bins;
  axisTitlePair[0]  = "#phi - #Psi_{2} (a.u.)"; 
  
  // delta eta
  const Int_t kNDeltaEtaBins = 80;
  Double_t deltaEtaBins[kNDeltaEtaBins+1];
  for(Int_t i = 0; i < kNDeltaEtaBins+1; i++)
    deltaEtaBins[i] = -2.0 + i * 0.05;
  iBinPair[1]       = kNDeltaEtaBins;
  dBinsPair[1]      = deltaEtaBins;
  axisTitlePair[1]  = "#Delta #eta"; 
  
  // delta phi
  const Int_t kNDeltaPhiBins = 72;
  Double_t deltaPhiBins[kNDeltaPhiBins+1];
  for(Int_t i = 0; i < kNDeltaPhiBins+1; i++){
    deltaPhiBins[i] = -180.0 + i * 5.;
  } 
  iBinPair[2]       = kNDeltaPhiBins;
  dBinsPair[2]      = deltaPhiBins;
  axisTitlePair[2]  = "#Delta #phi (#circ)"; 
  
  // pt(trigger-associated)
  const Int_t kNPtBins = 16;
  Double_t ptBins[kNPtBins+1] = {0.2,0.6,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,10.,12.,15.,20.};
  iBinSingle[1]       = kNPtBins;
  dBinsSingle[1]      = ptBins;
  axisTitleSingle[1]  = "p_{t}^{trig.} (GeV/c)"; 
  
  iBinPair[3]       = kNPtBins;
  dBinsPair[3]      = ptBins;
  axisTitlePair[3]  = "p_{t}^{trig.} (GeV/c)"; 
  
  iBinPair[4]       = kNPtBins;
  dBinsPair[4]      = ptBins;
  axisTitlePair[4]  = "p_{t}^{assoc.} (GeV/c)";   
  
  //=============================//
  //Generated: Single particle distributions
  fHistGeneratedPlus = new AliTHn("fHistGeneratedPlus",
				  "Generated positive primaries",
				  anaSteps,kVariablesSingle,iBinSingle);
  for (Int_t j = 0; j < kVariablesSingle; j++) {
    fHistGeneratedPlus->SetBinLimits(j, dBinsSingle[j]);
    fHistGeneratedPlus->SetVarTitle(j, axisTitleSingle[j]);
  }
  fOutputList->Add(fHistGeneratedPlus);

  fHistGeneratedMinus = new AliTHn("fHistGeneratedMinus",
				  "Generated positive primaries",
				  anaSteps,kVariablesSingle,iBinSingle);
  for (Int_t j = 0; j < kVariablesSingle; j++) {
    fHistGeneratedMinus->SetBinLimits(j, dBinsSingle[j]);
    fHistGeneratedMinus->SetVarTitle(j, axisTitleSingle[j]);
  }
  fOutputList->Add(fHistGeneratedMinus);

  //Survived: Single particle distributions
  fHistSurvivedPlus = new AliTHn("fHistSurvivedPlus",
				  "Survived positive primaries",
				  anaSteps,kVariablesSingle,iBinSingle);
  for (Int_t j = 0; j < kVariablesSingle; j++) {
    fHistSurvivedPlus->SetBinLimits(j, dBinsSingle[j]);
    fHistSurvivedPlus->SetVarTitle(j, axisTitleSingle[j]);
  }
  fOutputList->Add(fHistSurvivedPlus);

  fHistSurvivedMinus = new AliTHn("fHistSurvivedMinus",
				  "Survived positive primaries",
				  anaSteps,kVariablesSingle,iBinSingle);
  for (Int_t j = 0; j < kVariablesSingle; j++) {
    fHistSurvivedMinus->SetBinLimits(j, dBinsSingle[j]);
    fHistSurvivedMinus->SetVarTitle(j, axisTitleSingle[j]);
  }
  fOutputList->Add(fHistSurvivedMinus);

  //=============================//
  //Generated-Survived: Particle pairs +-
  fHistGeneratedPlusMinus = new AliTHn("fHistGeneratedPlusMinus",
				       "Generated +- primaries",
				       anaSteps,kVariablesPair,iBinPair);
  for (Int_t j = 0; j < kVariablesPair; j++) {
    fHistGeneratedPlusMinus->SetBinLimits(j, dBinsPair[j]);
    fHistGeneratedPlusMinus->SetVarTitle(j, axisTitlePair[j]);
  }
  fOutputList->Add(fHistGeneratedPlusMinus);

  fHistSurvivedPlusMinus = new AliTHn("fHistSurvivedPlusMinus",
				     "Survived +- primaries",
				     anaSteps,kVariablesPair,iBinPair);
  for (Int_t j = 0; j < kVariablesPair; j++) {
    fHistSurvivedPlusMinus->SetBinLimits(j, dBinsPair[j]);
    fHistSurvivedPlusMinus->SetVarTitle(j, axisTitlePair[j]);
  }
  fOutputList->Add(fHistSurvivedPlusMinus);

  //Generated-Survived: Particle pairs -+
  fHistGeneratedMinusPlus = new AliTHn("fHistGeneratedMinusPlus",
				       "Generated -+ primaries",
				       anaSteps,kVariablesPair,iBinPair);
  for (Int_t j = 0; j < kVariablesPair; j++) {
    fHistGeneratedMinusPlus->SetBinLimits(j, dBinsPair[j]);
    fHistGeneratedMinusPlus->SetVarTitle(j, axisTitlePair[j]);
  }
  fOutputList->Add(fHistGeneratedMinusPlus);
  
  fHistSurvivedMinusPlus = new AliTHn("fHistSurvivedMinusPlus",
				      "Survived -+ primaries",
				      anaSteps,kVariablesPair,iBinPair);
  for (Int_t j = 0; j < kVariablesPair; j++) {
    fHistSurvivedMinusPlus->SetBinLimits(j, dBinsPair[j]);
    fHistSurvivedMinusPlus->SetVarTitle(j, axisTitlePair[j]);
  }
  fOutputList->Add(fHistSurvivedMinusPlus);

  //Generated-Survived: Particle pairs ++
  fHistGeneratedPlusPlus = new AliTHn("fHistGeneratedPlusPlus",
				      "Generated ++ primaries",
				      anaSteps,kVariablesPair,iBinPair);
  for (Int_t j = 0; j < kVariablesPair; j++) {
    fHistGeneratedPlusPlus->SetBinLimits(j, dBinsPair[j]);
    fHistGeneratedPlusPlus->SetVarTitle(j, axisTitlePair[j]);
  }
  fOutputList->Add(fHistGeneratedPlusPlus);

  fHistSurvivedPlusPlus = new AliTHn("fHistSurvivedPlusPlus",
				     "Survived ++ primaries",
				     anaSteps,kVariablesPair,iBinPair);
  for (Int_t j = 0; j < kVariablesPair; j++) {
    fHistSurvivedPlusPlus->SetBinLimits(j, dBinsPair[j]);
    fHistSurvivedPlusPlus->SetVarTitle(j, axisTitlePair[j]);
  }
  fOutputList->Add(fHistSurvivedPlusPlus);

  //Generated-Survived: Particle pairs --
  fHistGeneratedMinusMinus = new AliTHn("fHistGeneratedMinusMinus",
					"Generated -- primaries",
					anaSteps,kVariablesPair,iBinPair);
  for (Int_t j = 0; j < kVariablesPair; j++) {
    fHistGeneratedMinusMinus->SetBinLimits(j, dBinsPair[j]);
    fHistGeneratedMinusMinus->SetVarTitle(j, axisTitlePair[j]);
  }
  fOutputList->Add(fHistGeneratedMinusMinus);
  
  fHistSurvivedMinusMinus = new AliTHn("fHistSurvivedMinusMinus",
				       "Survived -- primaries",
				       anaSteps,kVariablesPair,iBinPair);
  for (Int_t j = 0; j < kVariablesPair; j++) {
    fHistSurvivedMinusMinus->SetBinLimits(j, dBinsPair[j]);
    fHistSurvivedMinusMinus->SetVarTitle(j, axisTitlePair[j]);
  }
  fOutputList->Add(fHistSurvivedMinusMinus);
  //=============================//

  fQAList->Print();
  fOutputList->Print();

  PostData(1, fQAList);
  PostData(2, fOutputList);

  TH1::AddDirectory(oldStatus);

}

//________________________________________________________________________
void AliAnalysisTaskEfficiencyBFPsi::UserExec(Option_t *) {
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
    AliError("ERROR: Could not retrieve MC event");
    return;
  }
  AliStack* stack = mcEvent->Stack();
  if (!stack) {
    AliError("ERROR: Could not retrieve MC stack");
    return;
  }

  // arrays for 2 particle histograms
  Int_t nMCLabelCounter         = 0;
  const Int_t maxMCLabelCounter = 20000;

  Double_t phiMinusPsi[maxMCLabelCounter];
  Double_t eta[maxMCLabelCounter];
  Double_t pt[maxMCLabelCounter];
  Double_t phi[maxMCLabelCounter];
  Int_t level[maxMCLabelCounter];
  Int_t charge[maxMCLabelCounter];

  Double_t trackVariablesSingle[kVariablesSingle];
  Double_t trackVariablesPair[kVariablesPair];

  Double_t gReactionPlane       = 0.;
  AliGenHijingEventHeader* headerH = dynamic_cast<AliGenHijingEventHeader*>(dynamic_cast<AliMCEvent*>(mcEvent)->GenEventHeader());
  if (headerH) {
    gReactionPlane = headerH->ReactionPlaneAngle();
    gReactionPlane *= TMath::RadToDeg();
  }

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
	    if(TMath::Abs(vertex->GetX()) < fVxMax) {
	      if(TMath::Abs(vertex->GetY()) < fVyMax) {
		if(TMath::Abs(vertex->GetZ()) < fVzMax) {
		  fHistEventStats->Fill(4); //analyzed events
		  
		  Int_t nMCParticles = mcEvent->GetNumberOfTracks();
		  TArrayI labelMCArray(nMCParticles);
		  
		  for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) {
		    AliMCParticle *mcTrack = (AliMCParticle*) mcEvent->GetTrack(iTracks);
		    if (!mcTrack) {
		      AliError(Form("ERROR: Could not receive track %d (mc loop)", iTracks));
		      continue;
		    }
		    
		    //exclude particles generated out of the acceptance
		    Double_t vz = mcTrack->Zv();
		    if (TMath::Abs(vz) > 50.) continue;
		   
		    //acceptance
		    if(TMath::Abs(mcTrack->Eta()) > fEtaRangeMax) 
		      continue;
		    if((mcTrack->Pt() > fPtRangeMax)||(mcTrack->Pt() < fPtRangeMin)) 
		      continue;
		    //if((mcTrack->Phi() > fPhiRangeMax)||(mcTrack->Phi() < fPhiRangeMin)) 
		    //continue;
		    
		    TParticle* particle = mcTrack->Particle();
		    if(!particle) continue;
		    if(!stack->IsPhysicalPrimary(iTracks)) continue;

		    if(iTracks <= stack->GetNprimary()) {		      
		      Short_t gMCCharge = mcTrack->Charge();
		      Float_t firstPhi = mcTrack->Phi()*TMath::RadToDeg();
		      Float_t firstPt  = mcTrack->Pt();
    
		      // Event plane (determine psi bin)
		      Double_t gPsiMinusPhi    =   0.;
		      Double_t gPsiMinusPhiBin = -10.;
		      gPsiMinusPhi   = TMath::Abs(firstPhi - gReactionPlane);
		      //in-plane
		      if((gPsiMinusPhi <= 7.5)||
			 ((172.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 187.5)))
			gPsiMinusPhiBin = 0.0;
		      //intermediate
		      else if(((37.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 52.5))||
			      ((127.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 142.5))||
			      ((217.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 232.5))||
			      ((307.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 322.5)))
			gPsiMinusPhiBin = 1.0;
		      //out of plane
		      else if(((82.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 97.5))||
			      ((262.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 277.5)))
			gPsiMinusPhiBin = 2.0;
		      //everything else
		      else 
			gPsiMinusPhiBin = 3.0;

		      trackVariablesSingle[0]    =  gPsiMinusPhiBin;
		      trackVariablesSingle[1]    =  firstPt;  

		      if(gMCCharge > 0)
			fHistGeneratedPlus->Fill(trackVariablesSingle,0,1.);
		      else if(gMCCharge < 0)
			fHistGeneratedMinus->Fill(trackVariablesSingle,0,1.);
		      
		      labelMCArray.AddAt(iTracks,nMCLabelCounter);
		      if(nMCLabelCounter >= maxMCLabelCounter){
			AliWarning(Form("MC Label Counter > Limit (%d) --> stop loop here",maxMCLabelCounter));
			break;
		      }
		      
		      //fill the arrays for 2 particle analysis
		      phiMinusPsi[nMCLabelCounter]    = gPsiMinusPhiBin;
		      eta[nMCLabelCounter] = particle->Eta();
		      pt[nMCLabelCounter]  = particle->Pt();
		      phi[nMCLabelCounter] = particle->Phi()*TMath::RadToDeg();
		      charge[nMCLabelCounter] = gMCCharge;
		      // findable = generated in this case!
		      
		      level[nMCLabelCounter]  = 1;
		      nMCLabelCounter += 1;
		    }//primaries
		  }//loop over MC particles
		
		  fHistNMult->Fill(nMCLabelCounter);
		  
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
		  
		    Int_t mcGoods = nMCLabelCounter;
		    for (Int_t k = 0; k < mcGoods; k++) {
		      Int_t mcLabel = labelMCArray.At(k);
		      		    			      
		      if (mcLabel != TMath::Abs(label)) continue;
		      if(mcLabel != label) continue;
		      if(label > stack->GetNtrack()) continue;
		      
		      TParticle *particle = stack->Particle(label);
		      if(!particle) continue;
		      
		      //acceptance
		      if(TMath::Abs(particle->Eta()) > fEtaRangeMax) 
			continue;
		      if((particle->Pt() > fPtRangeMax)||(particle->Pt() <  fPtRangeMin)) 
			continue;
		      //if((particle->Phi() > fPhiRangeMax)||(particle->Phi() < fPhiRangeMin)) 
		      //continue;

		      if(!stack->IsPhysicalPrimary(label)) continue;
		      
		      if(label <= stack->GetNprimary()) {		
			Short_t gCharge = track->Charge();		
			// track cuts + analysis kinematic cuts
			if(fESDtrackCuts->AcceptTrack(track) && TMath::Abs(track->Eta()) < fMaxEta && track->Pt() > fMinPt && track->Pt() < fMaxPt ){
			  // survived
			  level[k]  = 2;

			  Float_t firstPhi = particle->Phi()*TMath::RadToDeg();
			  Float_t firstPt  = particle->Pt();
    
			  // Event plane (determine psi bin)
			  Double_t gPsiMinusPhi    =   0.;
			  Double_t gPsiMinusPhiBin = -10.;
			  gPsiMinusPhi = TMath::Abs(firstPhi - gReactionPlane);
			  //in-plane
			  if((gPsiMinusPhi <= 7.5)||
			     ((172.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 187.5)))
			    gPsiMinusPhiBin = 0.0;
			  //intermediate
			  else if(((37.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 52.5))||
				  ((127.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 142.5))||
				  ((217.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 232.5))||
				  ((307.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 322.5)))
			    gPsiMinusPhiBin = 1.0;
			  //out of plane
			  else if(((82.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 97.5))||
				  ((262.5 <= gPsiMinusPhi)&&(gPsiMinusPhi <= 277.5)))
			    gPsiMinusPhiBin = 2.0;
			  //everything else
			  else 
			    gPsiMinusPhiBin = 3.0;
			  
			  trackVariablesSingle[0]    =  gPsiMinusPhiBin;
			  trackVariablesSingle[1]    =  firstPt;  
			  
			  if(gCharge > 0)
			    fHistSurvivedPlus->Fill(trackVariablesSingle,0,1.);
			  else if(gCharge < 0)
			    fHistSurvivedMinus->Fill(trackVariablesSingle,0,1.);			  
			}//track cuts
		      }//primary particles
		    }//loop over the generated
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
    Float_t firstEta = eta[i];
    Float_t firstPhi = phi[i];
    Float_t firstPt  = pt[i];
    Float_t gPhisMinusPsiBin = phiMinusPsi[i];
    for (Int_t j = 0; j < nMCLabelCounter ; j++) {
      if(i == j) continue;

      Float_t secondEta = eta[j];
      Float_t secondPhi = phi[j];
      Float_t secondPt  = pt[j];

      trackVariablesPair[0]    =  gPhisMinusPsiBin;
      trackVariablesPair[1]    =  firstEta - secondEta;  // delta eta
      trackVariablesPair[2]    =  firstPhi - secondPhi;  // delta phi
      if (trackVariablesPair[2] > 180.)   // delta phi between -180 and 180 
	trackVariablesPair[2] -= 360.;
      if (trackVariablesPair[2] <  - 180.) 
	trackVariablesPair[2] += 360.;
      trackVariablesPair[3]    =  firstPt;      // pt trigger
      trackVariablesPair[4]    =  secondPt;  // pt

      //++ pairs
      if(charge[i] > 0 && charge[j] > 0 ) {
	if(level[i] > 0 && level[j] > 0) 
	  fHistGeneratedPlusPlus->Fill(trackVariablesPair,0,1.);
       
	if(level[i] > 1 && level[j] > 1) 
	  fHistSurvivedPlusPlus->Fill(trackVariablesPair,0,1.);
      }

      //-- pairs
      else if(charge[i] < 0 && charge[j] < 0 ) {
	if(level[i] > 0 && level[j] > 0) 
	  fHistGeneratedMinusMinus->Fill(trackVariablesPair,0,1.);
       
	if(level[i] > 1 && level[j] > 1) 
	  fHistSurvivedMinusMinus->Fill(trackVariablesPair,0,1.);
      }

      //+- pairs
      else if(charge[i] > 0 && charge[j] < 0 ) {
	if(level[i] > 0 && level[j] > 0) 
	  fHistGeneratedPlusMinus->Fill(trackVariablesPair,0,1.);
       
	if(level[i] > 1 && level[j] > 1) 
	  fHistSurvivedPlusMinus->Fill(trackVariablesPair,0,1.);
      }

      //-+ pairs
      else if(charge[i] < 0 && charge[j] > 0 ) {
	if(level[i] > 0 && level[j] > 0) 
	  fHistGeneratedMinusPlus->Fill(trackVariablesPair,0,1.);
       
	if(level[i] > 1 && level[j] > 1) 
	  fHistSurvivedMinusPlus->Fill(trackVariablesPair,0,1.);
      }
    }//second particle loop
  }//first particle loop
  
}

//________________________________________________________________________
void AliAnalysisTaskEfficiencyBFPsi::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}

//____________________________________________________________________//
Bool_t AliAnalysisTaskEfficiencyBFPsi::IsLabelUsed(TArrayI labelArray, Int_t label) {
  //Checks if the label is used already
  Bool_t status = kFALSE;
  for(Int_t i = 0; i < labelArray.GetSize(); i++) {
    if(labelArray.At(i) == label)
      status = kTRUE;
  }

  return status;
}
