/**************************************************************************
 * Author: Panos Christakoglou.                                           *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//-----------------------------------------------------------------
//           Balance Function class
//   This is the class to deal with the Balance Function analysis
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------


//ROOT
#include <Riostream.h>
#include <TMath.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TGraphErrors.h>
#include <TString.h>

#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"

#include "AliBalance.h"

ClassImp(AliBalance)

//____________________________________________________________________//
AliBalance::AliBalance() :
  TObject(), 
  bShuffle(kFALSE),
  fAnalysisLevel("ESD"),
  fAnalyzedEvents(0) ,
  fCentralityId(0) {
  // Default constructor
 
  for(Int_t i = 0; i < ANALYSIS_TYPES; i++){
    if(i == 6) {
      fNumberOfBins[i] = 180;
      fP1Start[i]      = -360.0;
      fP1Stop[i]       = 360.0;
      fP2Start[i]      = -360.0;
      fP2Stop[i]       = 360.0;
      fP2Step[i]       = 0.1;
    }
    else {
      fNumberOfBins[i] = 20;
      fP1Start[i]      = -1.0;
      fP1Stop[i]       = 1.0;
      fP2Start[i]      = 0.0;
      fP2Stop[i]       = 2.0;
    }
    fP2Step[i] = TMath::Abs(fP2Start - fP2Stop) / (Double_t)fNumberOfBins[i];
    fNn[i] = 0.0;
    fNp[i] = 0.0;

    for(Int_t j = 0; j < MAXIMUM_NUMBER_OF_STEPS; j++) {
      fNpp[i][j] = .0;
      fNnn[i][j] = .0;
      fNpn[i][j] = .0;
      fNnp[i][j] = .0;
      fB[i][j] = 0.0;
      ferror[i][j] = 0.0;
    }

    fHistP[i]  = NULL;
    fHistN[i]  = NULL;
    fHistPP[i] = NULL;
    fHistPN[i] = NULL;
    fHistNP[i] = NULL;
    fHistNN[i] = NULL;

  }
}


//____________________________________________________________________//
AliBalance::AliBalance(const AliBalance& balance):
  TObject(balance), bShuffle(balance.bShuffle), 
  fAnalysisLevel(balance.fAnalysisLevel),
  fAnalyzedEvents(balance.fAnalyzedEvents), 
  fCentralityId(balance.fCentralityId) {
  //copy constructor
  for(Int_t i = 0; i < ANALYSIS_TYPES; i++){
    fNn[i] = balance.fNn[i];
    fNp[i] = balance.fNp[i];

    fP1Start[i]      = balance.fP1Start[i];
    fP1Stop[i]       = balance.fP1Stop[i];
    fNumberOfBins[i] = balance.fNumberOfBins[i];
    fP2Start[i]      = balance.fP2Start[i];
    fP2Stop[i]       = balance.fP2Stop[i];
    fP2Step[i]       = balance.fP2Step[i];

    fHistP[i]        = balance.fHistP[i];
    fHistN[i]        = balance.fHistN[i];
    fHistPN[i]        = balance.fHistPN[i];
    fHistNP[i]        = balance.fHistNP[i];
    fHistPP[i]        = balance.fHistPP[i];
    fHistNN[i]        = balance.fHistNN[i];

    for(Int_t j = 0; j < MAXIMUM_NUMBER_OF_STEPS; j++) {
      fNpp[i][j] = .0;
      fNnn[i][j] = .0;
      fNpn[i][j] = .0;
      fNnp[i][j] = .0;
      fB[i][j] = 0.0;
      ferror[i][j] = 0.0;
    } 
  }
 }
 

//____________________________________________________________________//
AliBalance::~AliBalance() {
  // Destructor

  for(Int_t i = 0; i < ANALYSIS_TYPES; i++){
 
    delete fHistP[i];
    delete fHistN[i];
    delete fHistPN[i];
    delete fHistNP[i];
    delete fHistPP[i];
    delete fHistNN[i];
  
  }
}

//____________________________________________________________________//
void AliBalance::SetInterval(Int_t iAnalysisType,
			     Double_t p1Start, Double_t p1Stop,
			     Int_t ibins, Double_t p2Start, Double_t p2Stop) {
  // Sets the analyzed interval. 
  // Set the same Information for all analyses

  if(iAnalysisType == -1){             
    for(Int_t i = 0; i < ANALYSIS_TYPES; i++){
      fP1Start[i] = p1Start;
      fP1Stop[i] = p1Stop;
      fNumberOfBins[i] = ibins;
      fP2Start[i] = p2Start;
      fP2Stop[i] = p2Stop;
      fP2Step[i] = TMath::Abs(p2Start - p2Stop) / (Double_t)fNumberOfBins[i];
    }
  }
  // Set the Information for one analysis
  else if((iAnalysisType > -1) && (iAnalysisType < ANALYSIS_TYPES)) {
    fP1Start[iAnalysisType] = p1Start;
    fP1Stop[iAnalysisType] = p1Stop;
    fNumberOfBins[iAnalysisType] = ibins;
    fP2Start[iAnalysisType] = p2Start;
    fP2Stop[iAnalysisType] = p2Stop;
    fP2Step[iAnalysisType] = TMath::Abs(p2Start - p2Stop) / (Double_t)fNumberOfBins[iAnalysisType];
  }
  else {
    AliError("Wrong ANALYSIS number!");
  }
}

//____________________________________________________________________//
void AliBalance::InitHistograms() {
  //Initialize the histograms
  TString histName;
  for(Int_t iAnalysisType = 0; iAnalysisType < ANALYSIS_TYPES; iAnalysisType++) {
    histName = "fHistP"; histName += gBFAnalysisType[iAnalysisType]; 
    if(bShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistP[iAnalysisType] = new TH1D(histName.Data(),"",100,fP1Start[iAnalysisType],fP1Stop[iAnalysisType]);

    histName = "fHistN"; histName += gBFAnalysisType[iAnalysisType]; 
    if(bShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistN[iAnalysisType] = new TH1D(histName.Data(),"",100,fP1Start[iAnalysisType],fP1Stop[iAnalysisType]);
  
    histName = "fHistPN"; histName += gBFAnalysisType[iAnalysisType]; 
    if(bShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistPN[iAnalysisType] = new TH1D(histName.Data(),"",fNumberOfBins[iAnalysisType],fP2Start[iAnalysisType],fP2Stop[iAnalysisType]);
    
    histName = "fHistNP"; histName += gBFAnalysisType[iAnalysisType]; 
    if(bShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistNP[iAnalysisType] = new TH1D(histName.Data(),"",fNumberOfBins[iAnalysisType],fP2Start[iAnalysisType],fP2Stop[iAnalysisType]);
    
    histName = "fHistPP"; histName += gBFAnalysisType[iAnalysisType]; 
    if(bShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistPP[iAnalysisType] = new TH1D(histName.Data(),"",fNumberOfBins[iAnalysisType],fP2Start[iAnalysisType],fP2Stop[iAnalysisType]);
    
    histName = "fHistNN"; histName += gBFAnalysisType[iAnalysisType]; 
    if(bShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistNN[iAnalysisType] = new TH1D(histName.Data(),"",fNumberOfBins[iAnalysisType],fP2Start[iAnalysisType],fP2Stop[iAnalysisType]);
  }
}

//____________________________________________________________________//
void AliBalance::PrintAnalysisSettings() {
  
  Printf("======================================");
  Printf("Analysis level: %s",fAnalysisLevel.Data());
  Printf("======================================");
  for(Int_t ibin = 0; ibin < ANALYSIS_TYPES; ibin++){
    Printf("Interval info for variable %d",ibin);
    Printf("Analyzed interval (min.): %lf",fP2Start[ibin]);
    Printf("Analyzed interval (max.): %lf",fP2Stop[ibin]);
    Printf("Number of bins: %d",fNumberOfBins[ibin]);
    Printf("Step: %lf",fP2Step[ibin]);
    Printf("          ");
  }
  Printf("======================================");
}

//____________________________________________________________________//
void AliBalance::CalculateBalance(TObjArray *gTrackArray, vector<Int_t> chargeVector) {
  // Calculates the balance function
  fAnalyzedEvents++;
  Int_t i = 0 , j = 0;
  Int_t iBin = 0;
  
  AliVParticle* track = 0;
  AliVParticle* track1 = 0;
  AliVParticle* track2 = 0;

  // Initialize histograms if not done yet
  if(!fHistPN[0]){
    AliWarning("Histograms not yet initialized! --> Will be done now");
    AliWarning("This works only in local mode --> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
    InitHistograms();
  }

  //Printf("(AliBalance) Number of tracks: %d",gTrackArray->GetEntries());
  Int_t gNtrack = gTrackArray->GetEntries();
  for(i = 0; i < gNtrack; i++) {
    if(fAnalysisLevel == "ESD")
      track = dynamic_cast<AliESDtrack *>(gTrackArray->At(i));
    else if(fAnalysisLevel == "AOD")
      track = dynamic_cast<AliAODTrack *>(gTrackArray->At(i));
    else if(fAnalysisLevel == "MC")
      track = dynamic_cast<AliMCParticle *>(gTrackArray->At(i));

    if(track) {
      Short_t charge          = chargeVector.at(i);
      Double_t pseudorapidity = track->Eta();
      Double_t rapidity       = track->Y();
      Double_t phi            = track->Phi() * TMath::RadToDeg();
      
      //0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi
      for(Int_t iAnalysisType = 0; iAnalysisType < ANALYSIS_TYPES; iAnalysisType++) {
	if(iAnalysisType == kEta) {
	  if((pseudorapidity >= fP1Start[iAnalysisType]) && (pseudorapidity <= fP1Stop[iAnalysisType])) {
	    if(charge > 0) {
	      fNp[iAnalysisType] += 1.;
	      fHistP[iAnalysisType]->Fill(pseudorapidity);
	    }//charge > 0
	    if(charge < 0) {
	      fNn[iAnalysisType] += 1.;
	      fHistN[iAnalysisType]->Fill(pseudorapidity);
	    }//charge < 0
	  }//p1 interval check
	}//analysis type: eta
	else if(iAnalysisType == kPhi) {
	  if((phi >= fP1Start[iAnalysisType]) && (phi <= fP1Stop[iAnalysisType])) {
	    if(charge > 0) {
	      fNp[iAnalysisType] += 1.;
	      fHistP[iAnalysisType]->Fill(phi);
	    }//charge > 0
	    if(charge < 0) {
	      fNn[iAnalysisType] += 1.;
	      fHistN[iAnalysisType]->Fill(phi);
	    }//charge < 0
	  }//p1 interval check
	}//analysis type: phi
	else {
	  if((rapidity >= fP1Start[iAnalysisType]) && (rapidity <= fP1Stop[iAnalysisType])) {
	    if(charge > 0) {
	      fNp[iAnalysisType] += 1.;
	      fHistP[iAnalysisType]->Fill(rapidity);
	    }//charge > 0
	    if(charge < 0) {
	      fNn[iAnalysisType] += 1.;
	      fHistN[iAnalysisType]->Fill(rapidity);
	    }//charge < 0
	  }//p1 interval check
	}//analysis type: y, qside, qout, qlong, qinv
      }//analysis type loop
    }//track object valid
    else continue;
  }
  //Printf("Np: %lf - Nn: %lf",fNp[0],fNn[0]);

  Double_t dy = 0., deta = 0.;
  Double_t qLong = 0., qOut = 0., qSide = 0., qInv = 0.;
  Double_t dphi = 0.;

  Short_t charge1  = 0;
  Double_t p1      = 0.;
  Double_t pX1     = 0., pY1     = 0., pZ1     = 0.;
  Double_t eta1 = 0., rap1 = 0.;
  Double_t energy1 = 0.;
  Double_t phi1    = 0.;

  Short_t charge2  = 0;
  Double_t p2      = 0.;
  Double_t pX2     = 0., pY2     = 0., pZ2     = 0.;
  Double_t eta2 = 0., rap2 = 0.;
  Double_t energy2 = 0.;
  Double_t phi2    = 0.;
  //0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi
  for(i = 1; i < gNtrack; i++) {
    if(fAnalysisLevel == "ESD")
      track1 = dynamic_cast<AliESDtrack *>(gTrackArray->At(i));
    else if(fAnalysisLevel == "AOD")
      track1 = dynamic_cast<AliAODTrack *>(gTrackArray->At(i));
    else if(fAnalysisLevel == "MC")
      track1 = dynamic_cast<AliMCParticle *>(gTrackArray->At(i));


    if(track1) {
      charge1 = chargeVector.at(i);;
      p1      = track1->P();
      pX1     = track1->Px();
      pY1     = track1->Py();
      pZ1     = track1->Pz();
      energy1 = TMath::Sqrt(TMath::Power(track1->P(),2) +
			    TMath::Power(track1->M(),2));
      phi1    = track1->Phi() * TMath::RadToDeg();

    }
    else continue;
    
    for(j = 0; j < i; j++) {
      if(fAnalysisLevel == "ESD")
	track2 = dynamic_cast<AliESDtrack *>(gTrackArray->At(j));
      else if(fAnalysisLevel == "AOD")
	track2 = dynamic_cast<AliAODTrack *>(gTrackArray->At(j));
      else if(fAnalysisLevel == "MC")
	track2 = dynamic_cast<AliMCParticle *>(gTrackArray->At(j));
      
      if(track2) {
	charge2 = chargeVector.at(j);
	p2      = track2->P();
	pX2     = track2->Px();
	pY2     = track2->Py();
	pZ2     = track2->Pz();
	energy2 = TMath::Sqrt(TMath::Power(track2->P(),2) +
			      TMath::Power(track2->M(),2));
	phi2    = track2->Phi() * TMath::RadToDeg();

	// filling the arrays

	// RAPIDITY 
	rap1 = 0.5*log((energy1 + pZ1)/(energy1 - pZ1)); 
	rap2 = 0.5*log((energy2 + pZ2)/(energy2 - pZ2)); 
	dy = TMath::Abs(rap1 - rap2);

	// Eta
	eta1 = track1->Eta();
	eta2 = track2->Eta();
	deta = TMath::Abs(eta1 - eta2);

	//qlong
	Double_t eTot = energy1 + energy2;
	Double_t pxTot = pX1 + pX2;
	Double_t pyTot = pY1 + pY2;
	Double_t pzTot = pZ1 + pZ2;
	Double_t q0Tot = energy1 - energy2;
	Double_t qxTot = pX1 - pX2;
	Double_t qyTot = pY1 - pY2;
	Double_t qzTot = pZ1 - pZ2;
	Double_t snn = TMath::Power(eTot,2) - TMath::Power(pxTot,2) - TMath::Power(pyTot,2) - TMath::Power(pzTot,2);
	Double_t ptTot = TMath::Sqrt( TMath::Power(pxTot,2) + TMath::Power(pyTot,2));
	
	qLong = TMath::Abs(eTot*qzTot - pzTot*q0Tot)/TMath::Sqrt(snn + TMath::Power(ptTot,2));
	
	//qout
	qOut = TMath::Sqrt(snn/(snn + TMath::Power(ptTot,2))) * TMath::Abs(pxTot*qxTot + pyTot*qyTot)/ptTot;
	
	//qside
	qSide = TMath::Abs(pxTot*qyTot - pyTot*qxTot)/ptTot;
	
	//qinv
	qInv = TMath::Sqrt(TMath::Abs(-TMath::Power(q0Tot,2) +TMath::Power(qxTot,2) +TMath::Power(qyTot,2) +TMath::Power(qzTot,2)));
	
	//phi
	dphi = TMath::Abs(phi1 - phi2);

	//0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi
	for(Int_t iAnalysisType = 0; iAnalysisType < ANALYSIS_TYPES; iAnalysisType++) {
	  if(iAnalysisType == kRapidity) {
	    if((rap1 >= fP1Start[iAnalysisType]) && (rap1 <= fP1Stop[iAnalysisType]) && (rap2 >= fP1Start[iAnalysisType]) && (rap2 <= fP1Stop[iAnalysisType])) {
	      if( dy > fP2Start[kRapidity] && dy < fP2Stop[kRapidity]){
		iBin = Int_t((dy-fP2Start[kRapidity])/fP2Step[kRapidity]);
		if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
		  if((charge1 > 0)&&(charge2 > 0)) {
		    fNpp[kRapidity][iBin] += 1.;
		    fHistPP[kRapidity]->Fill(dy);
		  }
		  if((charge1 < 0)&&(charge2 < 0)) {
		    fNnn[kRapidity][iBin] += 1.;
		    fHistNN[kRapidity]->Fill(dy);
		  }
		  if((charge1 > 0)&&(charge2 < 0)) {
		    fNpn[kRapidity][iBin] += 1.;
		    fHistPN[kRapidity]->Fill(dy);
		  }
		  if((charge1 < 0)&&(charge2 > 0)) {
		    fNpn[kRapidity][iBin] += 1.;
		    fHistPN[kRapidity]->Fill(dy);
		  }
		}//BF binning check
	      }//p2 interval check
	    }//p1 interval check
	  }//rapidity
	  
	  if(iAnalysisType == kEta) {
	    if((eta1 >= fP1Start[iAnalysisType]) && (eta1 <= fP1Stop[iAnalysisType]) && (eta2 >= fP1Start[iAnalysisType]) && (eta2 <= fP1Stop[iAnalysisType])) {
	      if( deta > fP2Start[kEta] && deta < fP2Stop[kEta]){
		iBin = Int_t((deta-fP2Start[kEta])/fP2Step[kEta]);	
		if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
		  if((charge1 > 0)&&(charge2 > 0)) {
		    fNpp[kEta][iBin] += 1.;
		    fHistPP[kEta]->Fill(deta);
		  }
		  if((charge1 < 0)&&(charge2 < 0)) {
		    fNnn[kEta][iBin] += 1.;
		    fHistNN[kEta]->Fill(deta);
		  }
		  if((charge1 > 0)&&(charge2 < 0)) {
		    fNpn[kEta][iBin] += 1.;
		    fHistPN[kEta]->Fill(deta);
		  }
		  if((charge1 < 0)&&(charge2 > 0)) {
		    fNpn[kEta][iBin] += 1.;
		    fHistPN[kEta]->Fill(deta);
		  }
		}//BF binning check
	      }//p2 interval check
	    }//p1 interval check
	  }//pseudorapidity
	  
	  // Qlong, out, side, inv
	  // Check the p1 intervall for rapidity here (like for single tracks above)
	  if(iAnalysisType == kQlong) {
	    if((rap1 >= fP1Start[kRapidity]) && (rap1 <= fP1Stop[kRapidity]) && (rap2 >= fP1Start[kRapidity]) && (rap2 <= fP1Stop[kRapidity])) {
	      if( qLong > fP2Start[kQlong] && qLong < fP2Stop[kQlong]){
		iBin = Int_t((qLong-fP2Start[kQlong])/fP2Step[kQlong]);	
		if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
		  if((charge1 > 0)&&(charge2 > 0)) {
		    fNpp[kQlong][iBin] += 1.;
		    fHistPP[kQlong]->Fill(qLong);
		  }
		  if((charge1 < 0)&&(charge2 < 0)) {
		    fNnn[kQlong][iBin] += 1.;
		    fHistNN[kQlong]->Fill(qLong);
		  }
		  if((charge1 > 0)&&(charge2 < 0)) {
		    fNpn[kQlong][iBin] += 1.;
		    fHistPN[kQlong]->Fill(qLong);
		  }
		  if((charge1 < 0)&&(charge2 > 0)) {
		    fNpn[kQlong][iBin] += 1.;
		    fHistPN[kQlong]->Fill(qLong);
		  }
		}//BF binning check
	      }//p2 interval check
	    }//p1 interval check
	  }//qLong
	  
	  if(iAnalysisType == kQout) {
	    if((rap1 >= fP1Start[kRapidity]) && (rap1 <= fP1Stop[kRapidity]) && (rap2 >= fP1Start[kRapidity]) && (rap2 <= fP1Stop[kRapidity])) {
	      if( qOut > fP2Start[kQout] && qOut < fP2Stop[kQout]){
		iBin = Int_t((qOut-fP2Start[kQout])/fP2Step[kQout]);	
		if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
		  if((charge1 > 0)&&(charge2 > 0)) {
		    fNpp[kQout][iBin] += 1.;
		    fHistPP[kQout]->Fill(qOut);
		  }
		  if((charge1 < 0)&&(charge2 < 0)) {
		    fNnn[kQout][iBin] += 1.;
		    fHistNN[kQout]->Fill(qOut);
		  }
		  if((charge1 > 0)&&(charge2 < 0)) {
		    fNpn[kQout][iBin] += 1.;
		    fHistPN[kQout]->Fill(qOut);
		  }
		  if((charge1 < 0)&&(charge2 > 0)) {
		    fNpn[kQout][iBin] += 1.;
		    fHistPN[kQout]->Fill(qOut);
		  }
		}//BF binning check
	      }//p2 interval check
	    }//p1 interval check	
	  }//qOut
	
	  if(iAnalysisType == kQside) {
	    if((rap1 >= fP1Start[kRapidity]) && (rap1 <= fP1Stop[kRapidity]) && (rap2 >= fP1Start[kRapidity]) && (rap2 <= fP1Stop[kRapidity])) {
	      if( qSide > fP2Start[kQside] && qSide < fP2Stop[kQside]){
		iBin = Int_t((qSide-fP2Start[kQside])/fP2Step[kQside]);	
		if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
		  if((charge1 > 0)&&(charge2 > 0)) {
		    fNpp[kQside][iBin] += 1.;
		    fHistPP[kQside]->Fill(qSide);
		  }
		  if((charge1 < 0)&&(charge2 < 0)) {
		    fNnn[kQside][iBin] += 1.;
		    fHistNN[kQside]->Fill(qSide);
		  }
		  if((charge1 > 0)&&(charge2 < 0)) {
		    fNpn[kQside][iBin] += 1.;
		    fHistPN[kQside]->Fill(qSide);
		  }
		  if((charge1 < 0)&&(charge2 > 0)) {
		    fNpn[kQside][iBin] += 1.;
		    fHistPN[kQside]->Fill(qSide);
		  }
		}//BF binning check
	      }//p2 interval check
	    }//p1 interval check
	  }//qSide
	
	  if(iAnalysisType == kQinv) {
	    if((rap1 >= fP1Start[kRapidity]) && (rap1 <= fP1Stop[kRapidity]) && (rap2 >= fP1Start[kRapidity]) && (rap2 <= fP1Stop[kRapidity])) {
	      if( qInv > fP2Start[kQinv] && qInv < fP2Stop[kQinv]){
		iBin = Int_t((qInv-fP2Start[kQinv])/fP2Step[kQinv]);	
		if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
		  if((charge1 > 0)&&(charge2 > 0)) {
		    fNpp[kQinv][iBin] += 1.;
		  fHistPP[kQinv]->Fill(qInv);
		  }
		  if((charge1 < 0)&&(charge2 < 0)) {
		    fNnn[kQinv][iBin] += 1.;
		    fHistNN[kQinv]->Fill(qInv);
		  }
		  if((charge1 > 0)&&(charge2 < 0)) {
		    fNpn[kQinv][iBin] += 1.;
		    fHistPN[kQinv]->Fill(qInv);
		  }
		  if((charge1 < 0)&&(charge2 > 0)) {
		    fNpn[kQinv][iBin] += 1.;
		    fHistPN[kQinv]->Fill(qInv);
		  }
		}//BF binning check
	      }//p2 interval check
	    }//p1 interval check
	  }//qInv
	  
	  // Phi
	  if(iAnalysisType == kPhi) {
	    if((phi1 >= fP1Start[iAnalysisType]) && (phi1 <= fP1Stop[iAnalysisType]) && (phi2 >= fP1Start[iAnalysisType]) && (phi2 <= fP1Stop[iAnalysisType])) {
	      if( dphi > fP2Start[kPhi] && dphi < fP2Stop[kPhi]){
		iBin = Int_t((dphi-fP2Start[kPhi])/fP2Step[kPhi]);	
		if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
		  if((charge1 > 0)&&(charge2 > 0)) {
		    fNpp[kPhi][iBin] += 1.;
		    fHistPP[kPhi]->Fill(dphi);
		  }
		  if((charge1 < 0)&&(charge2 < 0)) {
		    fNnn[kPhi][iBin] += 1.;
		    fHistNN[kPhi]->Fill(dphi);
		  }
		  if((charge1 > 0)&&(charge2 < 0)) {
		    fNpn[kPhi][iBin] += 1.;
		    fHistPN[kPhi]->Fill(dphi);
		  }
		  if((charge1 < 0)&&(charge2 > 0)) {
		    fNpn[kPhi][iBin] += 1.;
		    fHistPN[kPhi]->Fill(dphi);
		  }
		}//BF binning check
	      }//p2 interval check
	    }//p1 interval check
	  }//phi
	}//analysis type loop
      }//track2 valid
      else continue;
    }//end of 2nd particle loop
  }//end of 1st particle loop
  //Printf("Number of analyzed events: %i",fAnalyzedEvents);
}  



//____________________________________________________________________//
Double_t AliBalance::GetBalance(Int_t iAnalysisType, Int_t p2) {
  // Returns the value of the balance function in bin p2
  fB[iAnalysisType][p2] = 0.5*(((fNpn[iAnalysisType][p2] - 2.*fNnn[iAnalysisType][p2])/fNn[iAnalysisType]) + ((fNpn[iAnalysisType][p2] - 2.*fNpp[iAnalysisType][p2])/fNp[iAnalysisType]))/fP2Step[iAnalysisType];
  
  return fB[iAnalysisType][p2];
}
    
//____________________________________________________________________//
Double_t AliBalance::GetError(Int_t iAnalysisType, Int_t p2) {        
  // Returns the error on the BF value for bin p2
  // The errors for fNn and fNp are neglected here (0.1 % of total error)
  /*ferror[iAnalysisType][p2] = TMath::Sqrt(Double_t(fNpp[iAnalysisType][p2])/(Double_t(fNp[iAnalysisType])*Double_t(fNp[iAnalysisType]))
			      + Double_t(fNnn[iAnalysisType][p2])/(Double_t(fNn[iAnalysisType])*Double_t(fNn[iAnalysisType]))
			      + Double_t(fNpn[iAnalysisType][p2])/(Double_t(fNp[iAnalysisType])*Double_t(fNp[iAnalysisType])) 
			      + Double_t(fNnp[iAnalysisType][p2])/(Double_t(fNp[iAnalysisType])*Double_t(fNp[iAnalysisType]))
			      //+ TMath::Power(fNpn[iAnalysisType][p2]-fNpp[iAnalysisType][p2],2)/TMath::Power(Double_t(fNp[iAnalysisType]),3)
			      //+ TMath::Power(fNnp[iAnalysisType][p2]-fNnn[iAnalysisType][p2],2)/TMath::Power(Double_t(fNn[iAnalysisType]),3) 
			       ) /fP2Step[iAnalysisType];*/

  ferror[iAnalysisType][p2] = TMath::Sqrt( Double_t(fNpp[iAnalysisType][p2])/(Double_t(fNp[iAnalysisType])*Double_t(fNp[iAnalysisType])) + 
					   Double_t(fNnn[iAnalysisType][p2])/(Double_t(fNn[iAnalysisType])*Double_t(fNn[iAnalysisType])) + 
					   Double_t(fNpn[iAnalysisType][p2])*TMath::Power((0.5/Double_t(fNp[iAnalysisType]) + 0.5/Double_t(fNn[iAnalysisType])),2))/fP2Step[iAnalysisType];
  
  return ferror[iAnalysisType][p2];
}
//____________________________________________________________________//
TGraphErrors *AliBalance::DrawBalance(Int_t iAnalysisType) {

  // Draws the BF
  Double_t x[MAXIMUM_NUMBER_OF_STEPS];
  Double_t xer[MAXIMUM_NUMBER_OF_STEPS];
  Double_t b[MAXIMUM_NUMBER_OF_STEPS];
  Double_t ber[MAXIMUM_NUMBER_OF_STEPS];

  if((fNp[iAnalysisType] == 0)||(fNn[iAnalysisType] == 0)) {
    cerr<<"Couldn't find any particles in the analyzed interval!!!"<<endl;
    return NULL;
  }
  
  for(Int_t i = 0; i < fNumberOfBins[iAnalysisType]; i++) {
    b[i] = GetBalance(iAnalysisType,i);
    ber[i] = GetError(iAnalysisType,i);
    x[i] = fP2Start[iAnalysisType] + fP2Step[iAnalysisType]*i + fP2Step[iAnalysisType]/2;
    xer[i] = 0.0;
  }
  
  TGraphErrors *gr = new TGraphErrors(fNumberOfBins[iAnalysisType],x,b,xer,ber);
  gr->GetXaxis()->SetTitleColor(1);
  if(iAnalysisType==0) {
    gr->SetTitle("Balance function B(#Delta y)");
    gr->GetXaxis()->SetTitle("#Delta y");
    gr->GetYaxis()->SetTitle("B(#Delta y)");
  }
  if(iAnalysisType==1) {
    gr->SetTitle("Balance function B(#Delta #eta)");
    gr->GetXaxis()->SetTitle("#Delta #eta");
    gr->GetYaxis()->SetTitle("B(#Delta #eta)");
  }
  if(iAnalysisType==2) {
    gr->SetTitle("Balance function B(q_{long})");
    gr->GetXaxis()->SetTitle("q_{long} (GeV/c)");
    gr->GetYaxis()->SetTitle("B(q_{long}) ((GeV/c)^{-1})");
  }
  if(iAnalysisType==3) {
    gr->SetTitle("Balance function B(q_{out})");
    gr->GetXaxis()->SetTitle("q_{out} (GeV/c)");
    gr->GetYaxis()->SetTitle("B(q_{out}) ((GeV/c)^{-1})");
  }
  if(iAnalysisType==4) {
    gr->SetTitle("Balance function B(q_{side})");
    gr->GetXaxis()->SetTitle("q_{side} (GeV/c)");
    gr->GetYaxis()->SetTitle("B(q_{side}) ((GeV/c)^{-1})");
  }
  if(iAnalysisType==5) {
    gr->SetTitle("Balance function B(q_{inv})");
    gr->GetXaxis()->SetTitle("q_{inv} (GeV/c)");
    gr->GetYaxis()->SetTitle("B(q_{inv}) ((GeV/c)^{-1})");
  }
  if(iAnalysisType==6) {
    gr->SetTitle("Balance function B(#Delta #phi)");
    gr->GetXaxis()->SetTitle("#Delta #phi");
    gr->GetYaxis()->SetTitle("B(#Delta #phi)");
  }

  return gr;
}

//____________________________________________________________________//
void AliBalance::PrintResults(Int_t iAnalysisType, TH1D *gHistBalance) {
  //Prints the calculated width of the BF and its error
  Double_t x[MAXIMUM_NUMBER_OF_STEPS];
  Double_t gSumXi = 0.0, gSumBi = 0.0, gSumBiXi = 0.0;
  Double_t gSumBiXi2 = 0.0, gSumBi2Xi2 = 0.0;
  Double_t gSumDeltaBi2 = 0.0, gSumXi2DeltaBi2 = 0.0;
  Double_t deltaBalP2 = 0.0, integral = 0.0;
  Double_t deltaErrorNew = 0.0;
  
  cout<<"=================================================="<<endl;
  for(Int_t i = 1; i <= fNumberOfBins[iAnalysisType]; i++) { 
    x[i-1] = fP2Start[iAnalysisType] + fP2Step[iAnalysisType]*i + fP2Step[iAnalysisType]/2;
    cout<<"B: "<<gHistBalance->GetBinContent(i)<<"\t Error: "<<gHistBalance->GetBinError(i)<<"\t bin: "<<gHistBalance->GetBinCenter(i)<<endl;
  } 
  cout<<"=================================================="<<endl;
  for(Int_t i = 2; i <= fNumberOfBins[iAnalysisType]; i++) {
    gSumXi += gHistBalance->GetBinCenter(i);
    gSumBi += gHistBalance->GetBinContent(i);
    gSumBiXi += gHistBalance->GetBinContent(i)*gHistBalance->GetBinCenter(i);
    gSumBiXi2 += gHistBalance->GetBinContent(i)*TMath::Power(gHistBalance->GetBinCenter(i),2);
    gSumBi2Xi2 += TMath::Power(gHistBalance->GetBinContent(i),2)*TMath::Power(gHistBalance->GetBinCenter(i),2);
    gSumDeltaBi2 +=  TMath::Power(gHistBalance->GetBinError(i),2);
    gSumXi2DeltaBi2 += TMath::Power(gHistBalance->GetBinCenter(i),2) * TMath::Power(gHistBalance->GetBinError(i),2);
    
    deltaBalP2 += fP2Step[iAnalysisType]*TMath::Power(gHistBalance->GetBinError(i),2);
    integral += fP2Step[iAnalysisType]*gHistBalance->GetBinContent(i);
  }
  for(Int_t i = 1; i < fNumberOfBins[iAnalysisType]; i++)
    deltaErrorNew += gHistBalance->GetBinError(i)*(gHistBalance->GetBinCenter(i)*gSumBi - gSumBiXi)/TMath::Power(gSumBi,2);
  
  Double_t integralError = TMath::Sqrt(deltaBalP2);
  
  Double_t delta = gSumBiXi / gSumBi;
  Double_t deltaError = (gSumBiXi / gSumBi) * TMath::Sqrt(TMath::Power((TMath::Sqrt(gSumXi2DeltaBi2)/gSumBiXi),2) + TMath::Power((gSumDeltaBi2/gSumBi),2) );
  
  cout<<"Width: "<<delta<<"\t Error: "<<deltaError<<endl;
  cout<<"New error: "<<deltaErrorNew<<endl;
  cout<<"Integral: "<<integral<<"\t Error: "<<integralError<<endl;
  cout<<"=================================================="<<endl;
}
 
//____________________________________________________________________//
TH1D *AliBalance::GetBalanceFunctionHistogram(Int_t iAnalysisType) {
  //Returns the BF histogram, extracted from the 6 TH1D objects 
  //(private members) of the AliBalance class.
  TString gAnalysisType[ANALYSIS_TYPES] = {"y","eta","qlong","qout","qside","qinv","phi"};
  TString histName = "gHistBalanceFunctionHistogram";
  histName += gAnalysisType[iAnalysisType];

  SetInterval(iAnalysisType, fHistP[iAnalysisType]->GetXaxis()->GetXmin(),
	      fHistP[iAnalysisType]->GetXaxis()->GetXmin(),
	      fHistPP[iAnalysisType]->GetNbinsX(),
	      fHistPP[iAnalysisType]->GetXaxis()->GetXmin(),
	      fHistPP[iAnalysisType]->GetXaxis()->GetXmax());

  TH1D *gHistBalanceFunctionHistogram = new TH1D(histName.Data(),"",fHistPP[iAnalysisType]->GetNbinsX(),fHistPP[iAnalysisType]->GetXaxis()->GetXmin(),fHistPP[iAnalysisType]->GetXaxis()->GetXmax());
  switch(iAnalysisType) {
  case kRapidity:
    gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("#Delta y");
    gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(#Delta y)");
    break;
  case kEta:
    gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("#Delta #eta");
    gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(#Delta #eta)");
    break;
  case kQlong:
    gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("q_{long} (GeV/c)");
    gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(q_{long})");
    break;
  case kQout:
    gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("q_{out} (GeV/c)");
    gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(q_{out})");
    break;
  case kQside:
    gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("q_{side} (GeV/c)");
    gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(q_{side})");
    break;
  case kQinv:
    gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("q_{inv} (GeV/c)");
    gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(q_{inv})");
    break;
  case kPhi:
    gHistBalanceFunctionHistogram->GetXaxis()->SetTitle("#Delta #phi (deg.)");
    gHistBalanceFunctionHistogram->GetYaxis()->SetTitle("B(#Delta #phi)");
    break;
  default:
    break;
  }
  
  TH1D *hTemp1 = dynamic_cast<TH1D *>(fHistPN[iAnalysisType]->Clone());
  TH1D *hTemp2 = dynamic_cast<TH1D *>(fHistPN[iAnalysisType]->Clone());
  TH1D *hTemp3 = dynamic_cast<TH1D *>(fHistNN[iAnalysisType]->Clone());
  TH1D *hTemp4 = dynamic_cast<TH1D *>(fHistPP[iAnalysisType]->Clone());

  if((hTemp1)&&(hTemp2)&&(hTemp3)&&(hTemp4)) {
    hTemp1->Sumw2();
    hTemp2->Sumw2();
    hTemp3->Sumw2();
    hTemp4->Sumw2();
    hTemp1->Add(hTemp3,-2.);
    hTemp1->Scale(1./fHistN[iAnalysisType]->GetEntries());
    hTemp2->Add(hTemp4,-2.);
    hTemp2->Scale(1./fHistP[iAnalysisType]->GetEntries());
    gHistBalanceFunctionHistogram->Add(hTemp1,hTemp2,1.,1.);
    gHistBalanceFunctionHistogram->Scale(0.5/fP2Step[iAnalysisType]);
  }
  
  PrintResults(iAnalysisType,gHistBalanceFunctionHistogram);

  return gHistBalanceFunctionHistogram;
}
