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

/* $Id: AliBalancePsi.cxx 54125 2012-01-24 21:07:41Z miweber $ */

//-----------------------------------------------------------------
//           Balance Function class
//   This is the class to deal with the Balance Function wrt Psi analysis
//   Origin: Panos Christakoglou, Nikhef, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------


//ROOT
#include <Riostream.h>
#include <TMath.h>
#include <TAxis.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TGraphErrors.h>
#include <TString.h>

#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"

#include "AliBalancePsi.h"

ClassImp(AliBalancePsi)

//____________________________________________________________________//
AliBalancePsi::AliBalancePsi() :
  TObject(), 
  bShuffle(kFALSE),
  fAnalysisLevel("ESD"),
  fAnalyzedEvents(0) ,
  fCentralityId(0) ,
  fCentStart(0.),
  fCentStop(0.),
  fPsiInterval(15.),
  fPsiNumberOfBins(24) {
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
    fCentStart = 0.;
    fCentStop  = 0.;

    for(Int_t k = 0; k < MAXIMUM_STEPS_IN_PSI; k++) {
      fNn[i][k] = 0.0;
      fNp[i][k] = 0.0;

      for(Int_t j = 0; j < MAXIMUM_NUMBER_OF_STEPS; j++) {
	fNpp[i][k][j] = .0;
	fNnn[i][k][j] = .0;
	fNpn[i][k][j] = .0;
	fNnp[i][k][j] = .0;
	fB[i][k][j] = 0.0;
	ferror[i][k][j] = 0.0;
      }
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
AliBalancePsi::AliBalancePsi(const AliBalancePsi& balance):
  TObject(balance), bShuffle(balance.bShuffle), 
  fAnalysisLevel(balance.fAnalysisLevel),
  fAnalyzedEvents(balance.fAnalyzedEvents), 
  fCentralityId(balance.fCentralityId),
  fCentStart(balance.fCentStart),
  fCentStop(balance.fCentStop),
  fPsiInterval(balance.fPsiInterval),
  fPsiNumberOfBins(balance.fPsiNumberOfBins) {
  //copy constructor
  for(Int_t i = 0; i < ANALYSIS_TYPES; i++){
    fP1Start[i]      = balance.fP1Start[i];
    fP1Stop[i]       = balance.fP1Stop[i];
    fNumberOfBins[i] = balance.fNumberOfBins[i];
    fP2Start[i]      = balance.fP2Start[i];
    fP2Stop[i]       = balance.fP2Stop[i];
    fP2Step[i]       = balance.fP2Step[i];
    fCentStart       = balance.fCentStart;
    fCentStop        = balance.fCentStop; 

    fHistP[i]        = balance.fHistP[i];
    fHistN[i]        = balance.fHistN[i];
    fHistPN[i]        = balance.fHistPN[i];
    fHistNP[i]        = balance.fHistNP[i];
    fHistPP[i]        = balance.fHistPP[i];
    fHistNN[i]        = balance.fHistNN[i];

    for(Int_t k = 0; k < MAXIMUM_STEPS_IN_PSI; k++) {
      fNn[i][k] = balance.fNn[i][k];
      fNp[i][k] = balance.fNp[i][k];
      
      for(Int_t j = 0; j < MAXIMUM_NUMBER_OF_STEPS; j++) {
        fNpp[i][k][j] = balance.fNpp[i][k][j];
	fNnn[i][k][j] = balance.fNnn[i][k][j];
	fNpn[i][k][j] = balance.fNpn[i][k][j];
	fNnp[i][k][j] = balance.fNnp[i][k][j];
	fB[i][k][j] = balance.fB[i][k][j];
	ferror[i][k][j] = balance.ferror[i][k][j];
      } 
    }
  }
}

//____________________________________________________________________//
AliBalancePsi::~AliBalancePsi() {
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
void AliBalancePsi::SetInterval(Int_t iAnalysisType,
				Double_t p1Start, Double_t p1Stop,
				Int_t ibins, 
				Double_t p2Start, Double_t p2Stop,
				Double_t psiInterval) {
  // Sets the analyzed interval. 
  // Set the same Information for all analyses
  fPsiInterval = psiInterval;
  fPsiNumberOfBins = (Int_t)(360./fPsiInterval);

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
void AliBalancePsi::InitHistograms() {
  //Initialize the histograms
  TString histName;
  for(Int_t iAnalysisType = 0; iAnalysisType < ANALYSIS_TYPES; iAnalysisType++) {
    histName = "fHistP"; histName += gBFPsiAnalysisType[iAnalysisType]; 
    if(bShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistP[iAnalysisType] = new TH3D(histName.Data(),"",(Int_t)(fCentStop-fCentStart),fCentStart,fCentStop,fPsiNumberOfBins,-fPsiInterval/2.,360.-fPsiInterval/2.,100,fP1Start[iAnalysisType],fP1Stop[iAnalysisType]);

    histName = "fHistN"; histName += gBFPsiAnalysisType[iAnalysisType]; 
    if(bShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistN[iAnalysisType] = new TH3D(histName.Data(),"",(Int_t)(fCentStop-fCentStart),fCentStart,fCentStop,fPsiNumberOfBins,-fPsiInterval/2.,360.-fPsiInterval/2.,100,fP1Start[iAnalysisType],fP1Stop[iAnalysisType]);
  
    histName = "fHistPN"; histName += gBFPsiAnalysisType[iAnalysisType]; 
    if(bShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistPN[iAnalysisType] = new TH3D(histName.Data(),"",(Int_t)(fCentStop-fCentStart),fCentStart,fCentStop,fPsiNumberOfBins,-fPsiInterval/2.,360.-fPsiInterval/2.,fNumberOfBins[iAnalysisType],fP2Start[iAnalysisType],fP2Stop[iAnalysisType]);
    
    histName = "fHistNP"; histName += gBFPsiAnalysisType[iAnalysisType]; 
    if(bShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistNP[iAnalysisType] = new TH3D(histName.Data(),"",(Int_t)(fCentStop-fCentStart),fCentStart,fCentStop,fPsiNumberOfBins,-fPsiInterval/2.,360.-fPsiInterval/2.,fNumberOfBins[iAnalysisType],fP2Start[iAnalysisType],fP2Stop[iAnalysisType]);
    
    histName = "fHistPP"; histName += gBFPsiAnalysisType[iAnalysisType]; 
    if(bShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistPP[iAnalysisType] = new TH3D(histName.Data(),"",(Int_t)(fCentStop-fCentStart),fCentStart,fCentStop,fPsiNumberOfBins,-fPsiInterval/2.,360.-fPsiInterval/2.,fNumberOfBins[iAnalysisType],fP2Start[iAnalysisType],fP2Stop[iAnalysisType]);
    
    histName = "fHistNN"; histName += gBFPsiAnalysisType[iAnalysisType]; 
    if(bShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistNN[iAnalysisType] = new TH3D(histName.Data(),"",(Int_t)(fCentStop-fCentStart),fCentStart,fCentStop,fPsiNumberOfBins,-fPsiInterval/2.,360.-fPsiInterval/2.,fNumberOfBins[iAnalysisType],fP2Start[iAnalysisType],fP2Stop[iAnalysisType]);
  }
}

//____________________________________________________________________//
void AliBalancePsi::PrintAnalysisSettings() {
  
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
void AliBalancePsi::CalculateBalance(Float_t fCentrality, 
				     Double_t gReactionPlane, 
				     vector<Double_t> **chargeVector) {
  // Calculates the balance function
  fAnalyzedEvents++;
  Int_t i = 0 , j = 0;
  Int_t iBin = 0;
  
  // Initialize histograms if not done yet
  if(!fHistPN[0]){
    AliWarning("Histograms not yet initialized! --> Will be done now");
    AliWarning("This works only in local mode --> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
    InitHistograms();
  }

  Int_t gNtrack = chargeVector[0]->size();
  //Printf("(AliBalancePsi) Number of tracks: %d",gNtrack);
  Double_t gPsiMinusPhi = 0.;
  Int_t gPsiBin = -10;
  for(i = 0; i < gNtrack;i++){
    Double_t charge          = chargeVector[0]->at(i);
    Double_t rapidity       = chargeVector[1]->at(i);
    Double_t pseudorapidity = chargeVector[2]->at(i);
    Double_t phi            = chargeVector[3]->at(i);
    gPsiMinusPhi   = phi - gReactionPlane;
    if(gPsiMinusPhi < -fPsiInterval/2) gPsiMinusPhi = 360. + gPsiMinusPhi;
    gPsiBin = (Int_t)((gPsiMinusPhi + 360./fPsiNumberOfBins/2.)/(360./fPsiNumberOfBins));
    //Printf("Psi: %lf - phi%lf - Delta phi: %lf - bin: %d",gReactionPlane,phi,gPsiMinusPhi,gPsiBin);
    //0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi
    for(Int_t iAnalysisType = 0; iAnalysisType < ANALYSIS_TYPES; iAnalysisType++) {
      if(iAnalysisType == kEta) {
	if((pseudorapidity >= fP1Start[iAnalysisType]) && (pseudorapidity <= fP1Stop[iAnalysisType])) {
	  if(charge > 0) {
	    fNp[iAnalysisType][gPsiBin] += 1.;
	    fHistP[iAnalysisType]->Fill(fCentrality,gPsiMinusPhi,pseudorapidity);
	  }//charge > 0
	  if(charge < 0) {
	    fNn[iAnalysisType][gPsiBin] += 1.;
	    fHistN[iAnalysisType]->Fill(fCentrality,gPsiMinusPhi,pseudorapidity);
	  }//charge < 0
	}//p1 interval check
      }//analysis type: eta
      else if(iAnalysisType == kPhi) {
	if((phi >= fP1Start[iAnalysisType]) && (phi <= fP1Stop[iAnalysisType])) {
	  if(charge > 0) {
	    fNp[iAnalysisType][gPsiBin] += 1.;
	    fHistP[iAnalysisType]->Fill(fCentrality,gPsiMinusPhi,phi);
	  }//charge > 0
	  if(charge < 0) {
	    fNn[iAnalysisType][gPsiBin] += 1.;
	    fHistN[iAnalysisType]->Fill(fCentrality,gPsiMinusPhi,phi);
	  }//charge < 0
	}//p1 interval check
      }//analysis type: phi
      else {
	if((rapidity >= fP1Start[iAnalysisType]) && (rapidity <= fP1Stop[iAnalysisType])) {
	  if(charge > 0) {
	    fNp[iAnalysisType][gPsiBin] += 1.;
	    fHistP[iAnalysisType]->Fill(fCentrality,gPsiMinusPhi,rapidity);
	  }//charge > 0
	  if(charge < 0) {
	    fNn[iAnalysisType][gPsiBin] += 1.;
	    fHistN[iAnalysisType]->Fill(fCentrality,gPsiMinusPhi,rapidity);
	  }//charge < 0
	}//p1 interval check
      }//analysis type: y, qside, qout, qlong, qinv
    }//analysis type loop
  }
  
  //Printf("Np: %lf - Nn: %lf",fNp[0],fNn[0]);
  
  Double_t dy = 0., deta = 0.;
  Double_t qLong = 0., qOut = 0., qSide = 0., qInv = 0.;
  Double_t dphi = 0.;

  Double_t charge1  = 0;
  Double_t eta1 = 0., rap1 = 0.;
  Double_t px1 = 0., py1 = 0., pz1 = 0.;
  Double_t pt1 = 0.;
  Double_t energy1 = 0.;
  Double_t phi1    = 0.;

  Double_t charge2  = 0;
  Double_t eta2 = 0., rap2 = 0.;
  Double_t px2 = 0., py2 = 0., pz2 = 0.;
  Double_t pt2 = 0.;
  Double_t energy2 = 0.;
  Double_t phi2    = 0.;
  //0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi
  for(i = 1; i < gNtrack; i++) {

      charge1 = chargeVector[0]->at(i);
      rap1    = chargeVector[1]->at(i);
      eta1    = chargeVector[2]->at(i);
      phi1    = chargeVector[3]->at(i);
      px1     = chargeVector[4]->at(i);
      py1     = chargeVector[5]->at(i);
      pz1     = chargeVector[6]->at(i);
      pt1     = chargeVector[7]->at(i);
      energy1 = chargeVector[8]->at(i);
      gPsiMinusPhi   = phi1 - gReactionPlane;
      //if(phi1 > 180.) phi1 -= 360.;
      if(gPsiMinusPhi < -fPsiInterval/2) gPsiMinusPhi = 360. + gPsiMinusPhi;
      gPsiBin = (Int_t)((gPsiMinusPhi + 360./fPsiNumberOfBins/2.)/(360./fPsiNumberOfBins));

      for(j = 0; j < i; j++) {
	charge2 = chargeVector[0]->at(j);
	rap2    = chargeVector[1]->at(j);
	eta2    = chargeVector[2]->at(j);
	phi2    = chargeVector[3]->at(j);
	px2     = chargeVector[4]->at(j);
	py2     = chargeVector[5]->at(j);
	pz2     = chargeVector[6]->at(j);
	pt2     = chargeVector[7]->at(i);
	energy2 = chargeVector[8]->at(j);
	//if(phi2 > 180.) phi2 -= 360.;
	
	// filling the arrays
	// RAPIDITY 
	dy = rap1 - rap2;
	
	// Eta
	deta = eta1 - eta2;
	
	//qlong
	Double_t eTot = energy1 + energy2;
	Double_t pxTot = px1 + px2;
	Double_t pyTot = py1 + py2;
	Double_t pzTot = pz1 + pz2;
	Double_t q0Tot = energy1 - energy2;
	Double_t qxTot = px1 - px2;
	Double_t qyTot = py1 - py2;
	Double_t qzTot = pz1 - pz2;
	
	Double_t eTot2 = eTot*eTot;
	Double_t pTot2 = pxTot*pxTot + pyTot*pyTot + pzTot*pzTot;
	Double_t pzTot2 = pzTot*pzTot;

	Double_t q0Tot2 = q0Tot*q0Tot;
	Double_t qTot2  = qxTot*qxTot + qyTot*qyTot + qzTot*qzTot;

	Double_t snn    = eTot2 - pTot2;
	Double_t ptTot2 = pTot2 - pzTot2 ;
	Double_t ptTot  = TMath::Sqrt( ptTot2 );
	
	qLong = TMath::Abs(eTot*qzTot - pzTot*q0Tot)/TMath::Sqrt(snn + ptTot2);
	
	//qout
	qOut = TMath::Sqrt(snn/(snn + ptTot2)) * TMath::Abs(pxTot*qxTot + pyTot*qyTot)/ptTot;
	
	//qside
	qSide = TMath::Abs(pxTot*qyTot - pyTot*qxTot)/ptTot;
	
	//qinv
	qInv = TMath::Sqrt(TMath::Abs(-q0Tot2 + qTot2 ));
	
	//phi
	dphi = phi2 - phi1;
	if(dphi < -180.) dphi = -360. - dphi;  //dphi should be between -180 and 180!
	else if(dphi > 180.) dphi = 360. - dphi;  //dphi should be between -180 and 180!

	//0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi
	if((rap1 >= fP1Start[kRapidity]) && (rap1 <= fP1Stop[kRapidity]) && (rap2 >= fP1Start[kRapidity]) && (rap2 <= fP1Stop[kRapidity])) {

	  // rapidity
	  if( dy > fP2Start[kRapidity] && dy < fP2Stop[kRapidity]){
	    iBin = Int_t((dy-fP2Start[kRapidity])/fP2Step[kRapidity]);
	    if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
	      
	      if((charge1 > 0)&&(charge2 > 0)) {
		fNpp[kRapidity][gPsiBin][iBin] += 1.;
		fHistPP[kRapidity]->Fill(fCentrality,gPsiMinusPhi,dy);
	      }
	      else if((charge1 < 0)&&(charge2 < 0)) {
		fNnn[kRapidity][gPsiBin][iBin] += 1.;
		fHistNN[kRapidity]->Fill(fCentrality,gPsiMinusPhi,dy);
	      }
	      else if((charge1 > 0)&&(charge2 < 0)) {
		fNpn[kRapidity][gPsiBin][iBin] += 1.;
		fHistPN[kRapidity]->Fill(fCentrality,gPsiMinusPhi,dy);
	      }
	      else if((charge1 < 0)&&(charge2 > 0)) {
		fNpn[kRapidity][gPsiBin][iBin] += 1.;
		fHistPN[kRapidity]->Fill(fCentrality,gPsiMinusPhi,dy);
	      }
	    }//BF binning check
	  }//p2 interval check
	}//p1 interval check
	
	// pseudorapidity
	if((eta1 >= fP1Start[kEta]) && (eta1 <= fP1Stop[kEta]) && (eta2 >= fP1Start[kEta]) && (eta2 <= fP1Stop[kEta])) {
	  if( deta > fP2Start[kEta] && deta < fP2Stop[kEta]){
	    iBin = Int_t((deta-fP2Start[kEta])/fP2Step[kEta]);	
	    if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
	      if((charge1 > 0)&&(charge2 > 0)) {
		fNpp[kEta][gPsiBin][iBin] += 1.;
		fHistPP[kEta]->Fill(fCentrality,gPsiMinusPhi,deta);
	      }
	      if((charge1 < 0)&&(charge2 < 0)) {
		fNnn[kEta][gPsiBin][iBin] += 1.;
		fHistNN[kEta]->Fill(fCentrality,gPsiMinusPhi,deta);
	      }
	      if((charge1 > 0)&&(charge2 < 0)) {
		fNpn[kEta][gPsiBin][iBin] += 1.;
		fHistPN[kEta]->Fill(fCentrality,gPsiMinusPhi,deta);
	      }
	      if((charge1 < 0)&&(charge2 > 0)) {
		fNpn[kEta][gPsiBin][iBin] += 1.;
		fHistPN[kEta]->Fill(fCentrality,gPsiMinusPhi,deta);
	      }
	    }//BF binning check
	  }//p2 interval check
	}//p1 interval check
	
	// Qlong, out, side, inv
	// Check the p1 intervall for rapidity here (like for single tracks above)
	if((rap1 >= fP1Start[kRapidity]) && (rap1 <= fP1Stop[kRapidity]) && (rap2 >= fP1Start[kRapidity]) && (rap2 <= fP1Stop[kRapidity])) {
	  if( qLong > fP2Start[kQlong] && qLong < fP2Stop[kQlong]){
	    iBin = Int_t((qLong-fP2Start[kQlong])/fP2Step[kQlong]);	
	    if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
	      if((charge1 > 0)&&(charge2 > 0)) {
		fNpp[kQlong][gPsiBin][iBin] += 1.;
		fHistPP[kQlong]->Fill(fCentrality,gPsiMinusPhi,qLong);
	      }
	      if((charge1 < 0)&&(charge2 < 0)) {
		fNnn[kQlong][gPsiBin][iBin] += 1.;
		fHistNN[kQlong]->Fill(fCentrality,gPsiMinusPhi,qLong);
	      }
	      if((charge1 > 0)&&(charge2 < 0)) {
		fNpn[kQlong][gPsiBin][iBin] += 1.;
		fHistPN[kQlong]->Fill(fCentrality,gPsiMinusPhi,qLong);
	      }
	      if((charge1 < 0)&&(charge2 > 0)) {
		fNpn[kQlong][gPsiBin][iBin] += 1.;
		fHistPN[kQlong]->Fill(fCentrality,gPsiMinusPhi,qLong);
	      }
	    }//BF binning check
	  }//p2 interval check
	}//p1 interval check
	  
	if((rap1 >= fP1Start[kRapidity]) && (rap1 <= fP1Stop[kRapidity]) && (rap2 >= fP1Start[kRapidity]) && (rap2 <= fP1Stop[kRapidity])) {
	  if( qOut > fP2Start[kQout] && qOut < fP2Stop[kQout]){
	    iBin = Int_t((qOut-fP2Start[kQout])/fP2Step[kQout]);	
	    if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
	      if((charge1 > 0)&&(charge2 > 0)) {
		fNpp[kQout][gPsiBin][iBin] += 1.;
		fHistPP[kQout]->Fill(fCentrality,gPsiMinusPhi,qOut);
	  	  }
	      if((charge1 < 0)&&(charge2 < 0)) {
		fNnn[kQout][gPsiBin][iBin] += 1.;
		fHistNN[kQout]->Fill(fCentrality,gPsiMinusPhi,qOut);
	      }
	      if((charge1 > 0)&&(charge2 < 0)) {
		fNpn[kQout][gPsiBin][iBin] += 1.;
		fHistPN[kQout]->Fill(fCentrality,gPsiMinusPhi,qOut);
	      }
	      if((charge1 < 0)&&(charge2 > 0)) {
		fNpn[kQout][gPsiBin][iBin] += 1.;
		fHistPN[kQout]->Fill(fCentrality,gPsiMinusPhi,qOut);
	      }
	    }//BF binning check
	  }//p2 interval check
	}//p1 interval check	
	
	if((rap1 >= fP1Start[kRapidity]) && (rap1 <= fP1Stop[kRapidity]) && (rap2 >= fP1Start[kRapidity]) && (rap2 <= fP1Stop[kRapidity])) {
	  if( qSide > fP2Start[kQside] && qSide < fP2Stop[kQside]){
	    iBin = Int_t((qSide-fP2Start[kQside])/fP2Step[kQside]);	
	    if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
	      if((charge1 > 0)&&(charge2 > 0)) {
		fNpp[kQside][gPsiBin][iBin] += 1.;
		fHistPP[kQside]->Fill(fCentrality,gPsiMinusPhi,qSide);
	      }
	      if((charge1 < 0)&&(charge2 < 0)) {
		fNnn[kQside][gPsiBin][iBin] += 1.;
		fHistNN[kQside]->Fill(fCentrality,gPsiMinusPhi,qSide);
	      }
	      if((charge1 > 0)&&(charge2 < 0)) {
		fNpn[kQside][gPsiBin][iBin] += 1.;
		fHistPN[kQside]->Fill(fCentrality,gPsiMinusPhi,qSide);
	      }
	      if((charge1 < 0)&&(charge2 > 0)) {
		fNpn[kQside][gPsiBin][iBin] += 1.;
		fHistPN[kQside]->Fill(fCentrality,gPsiMinusPhi,qSide);
		  	  }
	    }//BF binning check
	  }//p2 interval check
	}//p1 interval check
	
	if((rap1 >= fP1Start[kRapidity]) && (rap1 <= fP1Stop[kRapidity]) && (rap2 >= fP1Start[kRapidity]) && (rap2 <= fP1Stop[kRapidity])) {
	  if( qInv > fP2Start[kQinv] && qInv < fP2Stop[kQinv]){
	    iBin = Int_t((qInv-fP2Start[kQinv])/fP2Step[kQinv]);	
	    if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
	      if((charge1 > 0)&&(charge2 > 0)) {
		fNpp[kQinv][gPsiBin][iBin] += 1.;
		fHistPP[kQinv]->Fill(fCentrality,gPsiMinusPhi,qInv);
	      }
	      if((charge1 < 0)&&(charge2 < 0)) {
		fNnn[kQinv][gPsiBin][iBin] += 1.;
		fHistNN[kQinv]->Fill(fCentrality,gPsiMinusPhi,qInv);
	      }
	      if((charge1 > 0)&&(charge2 < 0)) {
		fNpn[kQinv][gPsiBin][iBin] += 1.;
		fHistPN[kQinv]->Fill(fCentrality,gPsiMinusPhi,qInv);
	      }
	      if((charge1 < 0)&&(charge2 > 0)) {
		fNpn[kQinv][gPsiBin][iBin] += 1.;
		fHistPN[kQinv]->Fill(fCentrality,gPsiMinusPhi,qInv);
	      }
	    }//BF binning check
	  }//p2 interval check
	}//p1 interval check
	
	// Phi
	if((phi1 >= fP1Start[kPhi]) && (phi1 <= fP1Stop[kPhi]) && (phi2 >= fP1Start[kPhi]) && (phi2 <= fP1Stop[kPhi])) {
	  if( dphi > fP2Start[kPhi] && dphi < fP2Stop[kPhi]){
	    iBin = Int_t((dphi-fP2Start[kPhi])/fP2Step[kPhi]);	
	    if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
	      if((charge1 > 0)&&(charge2 > 0)) {
		fNpp[kPhi][gPsiBin][iBin] += 1.;
		fHistPP[kPhi]->Fill(fCentrality,gPsiMinusPhi,dphi);
	      }
	      if((charge1 < 0)&&(charge2 < 0)) {
		fNnn[kPhi][gPsiBin][iBin] += 1.;
		fHistNN[kPhi]->Fill(fCentrality,gPsiMinusPhi,dphi);
	      }
	      if((charge1 > 0)&&(charge2 < 0)) {
		fNpn[kPhi][gPsiBin][iBin] += 1.;
		fHistPN[kPhi]->Fill(fCentrality,gPsiMinusPhi,dphi);
	      }
	      if((charge1 < 0)&&(charge2 > 0)) {
		fNpn[kPhi][gPsiBin][iBin] += 1.;
		fHistPN[kPhi]->Fill(fCentrality,gPsiMinusPhi,dphi);
	      }
	    }//BF binning check
	  }//p2 interval check
	}//p1 interval check
    }//end of 2nd particle loop
  }//end of 1st particle loop
  //Printf("Number of analyzed events: %i",fAnalyzedEvents);
  //Printf("DeltaEta NN[0] = %.0f, PP[0] = %.0f, NP[0] = %.0f, PN[0] = %.0f",fNnn[kEta][0],fNpp[kEta][0],fNnp[kEta][0],fNpn[kEta][0]);
}  


//____________________________________________________________________//
Double_t AliBalancePsi::GetBalance(Int_t iAnalysisType, Int_t gPsiBin, 
				   Int_t p2) {
  // Returns the value of the balance function in bin p2
  fB[iAnalysisType][gPsiBin][p2] = 0.5*(((fNpn[iAnalysisType][gPsiBin][p2] - 2.*fNnn[iAnalysisType][gPsiBin][p2])/fNn[iAnalysisType][gPsiBin]) + ((fNpn[iAnalysisType][gPsiBin][p2] - 2.*fNpp[iAnalysisType][gPsiBin][p2])/fNp[iAnalysisType][gPsiBin]))/fP2Step[iAnalysisType];
  
  return fB[iAnalysisType][gPsiBin][p2];
}
    
//____________________________________________________________________//
Double_t AliBalancePsi::GetError(Int_t iAnalysisType, Int_t gPsiBin, 
				 Int_t p2) {        
  // Returns the error on the BF value for bin p2
  // The errors for fNn and fNp are neglected here (0.1 % of total error)
  ferror[iAnalysisType][gPsiBin][p2] = TMath::Sqrt( Double_t(fNpp[iAnalysisType][gPsiBin][p2])/(Double_t(fNp[iAnalysisType][gPsiBin])*Double_t(fNp[iAnalysisType][gPsiBin])) + 
					   Double_t(fNnn[iAnalysisType][gPsiBin][p2])/(Double_t(fNn[iAnalysisType][gPsiBin])*Double_t(fNn[iAnalysisType][gPsiBin])) + 
					   Double_t(fNpn[iAnalysisType][gPsiBin][p2])*TMath::Power((0.5/Double_t(fNp[iAnalysisType][gPsiBin]) + 0.5/Double_t(fNn[iAnalysisType][gPsiBin])),2))/fP2Step[iAnalysisType];
  
  return ferror[iAnalysisType][gPsiBin][p2];
}
//____________________________________________________________________//
TGraphErrors *AliBalancePsi::DrawBalance(Int_t iAnalysisType, Int_t gPsiBin) {
  // Draws the BF
  Double_t x[MAXIMUM_NUMBER_OF_STEPS];
  Double_t xer[MAXIMUM_NUMBER_OF_STEPS];
  Double_t b[MAXIMUM_NUMBER_OF_STEPS];
  Double_t ber[MAXIMUM_NUMBER_OF_STEPS];

  if((fNp[iAnalysisType][gPsiBin] == 0)||(fNn[iAnalysisType][gPsiBin] == 0)) {
    cerr<<"Couldn't find any particles in the analyzed interval!!!"<<endl;
    return NULL;
  }
  
  for(Int_t i = 0; i < fNumberOfBins[iAnalysisType]; i++) {
    b[i] = GetBalance(iAnalysisType,gPsiBin,i);
    ber[i] = GetError(iAnalysisType,gPsiBin,i);
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
void AliBalancePsi::PrintResults(Int_t iAnalysisType, TH1D *gHistBalance) {
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
    //cout<<"B: "<<gHistBalance->GetBinContent(i)<<"\t Error: "<<gHistBalance->GetBinError(i)<<"\t bin: "<<gHistBalance->GetBinCenter(i)<<endl;
  } 
  //cout<<"=================================================="<<endl;
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
  cout<<"Analysis type: "<<gBFPsiAnalysisType[iAnalysisType].Data()<<endl;
  cout<<"Width: "<<delta<<"\t Error: "<<deltaError<<endl;
  cout<<"New error: "<<deltaErrorNew<<endl;
  cout<<"Integral: "<<integral<<"\t Error: "<<integralError<<endl;
  cout<<"=================================================="<<endl;
}
 
//____________________________________________________________________//
TH1D *AliBalancePsi::GetBalanceFunctionHistogram(Int_t iAnalysisType,Double_t centrMin, Double_t centrMax, Double_t psiMin, Double_t psiMax) {
  //Returns the BF histogram, extracted from the 6 TH2D objects 
  //(private members) of the AliBalancePsi class.
  //
  TString gAnalysisType[ANALYSIS_TYPES] = {"y","eta","qlong","qout","qside","qinv","phi"};
  TString histName = "gHistBalanceFunctionHistogram";
  histName += gAnalysisType[iAnalysisType];

  SetInterval(iAnalysisType, fHistP[iAnalysisType]->GetZaxis()->GetXmin(),
	      fHistP[iAnalysisType]->GetZaxis()->GetXmin(),
	      fHistPP[iAnalysisType]->GetNbinsZ(),
	      fHistPP[iAnalysisType]->GetZaxis()->GetXmin(),
	      fHistPP[iAnalysisType]->GetZaxis()->GetXmax(),
	      psiMax-psiMin);

  // determine the projection thresholds
  Int_t binMinX, binMinY, binMinZ;
  Int_t binMaxX, binMaxY, binMaxZ;

  fHistPP[iAnalysisType]->GetBinXYZ(fHistPP[iAnalysisType]->FindBin(centrMin,psiMin),binMinX,binMinY,binMinZ);
  fHistPP[iAnalysisType]->GetBinXYZ(fHistPP[iAnalysisType]->FindBin(centrMax,psiMax),binMaxX,binMaxY,binMaxZ);

  TH1D *gHistBalanceFunctionHistogram = new TH1D(histName.Data(),"",fHistPP[iAnalysisType]->GetNbinsZ(),fHistPP[iAnalysisType]->GetZaxis()->GetXmin(),fHistPP[iAnalysisType]->GetZaxis()->GetXmax());
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

  TH1D *hTemp1 = dynamic_cast<TH1D *>(fHistPN[iAnalysisType]->ProjectionZ(Form("%s_Cent_%.0f_%.0f_Psi_%.0f_%.0f",fHistPN[iAnalysisType]->GetName(),centrMin,centrMax,psiMin,psiMax),binMinX,binMaxX,binMinY,binMaxY));
  TH1D *hTemp2 = dynamic_cast<TH1D *>(fHistPN[iAnalysisType]->ProjectionZ(Form("%s_Cent_%.0f_%.0f_Psi_%.0f_%.0f_copy",fHistPN[iAnalysisType]->GetName(),centrMin,centrMax,psiMin,psiMax),binMinX,binMaxX,binMinY,binMaxY));
  TH1D *hTemp3 = dynamic_cast<TH1D *>(fHistNN[iAnalysisType]->ProjectionZ(Form("%s_Cent_%.0f_%.0f_Psi_%.0f_%.0f",fHistNN[iAnalysisType]->GetName(),centrMin,centrMax,psiMin,psiMax),binMinX,binMaxX,binMinY,binMaxY));
  TH1D *hTemp4 = dynamic_cast<TH1D *>(fHistPP[iAnalysisType]->ProjectionZ(Form("%s_Cent_%.0f_%.0f_Psi_%.0f_%.0f",fHistPP[iAnalysisType]->GetName(),centrMin,centrMax,psiMin,psiMax),binMinX,binMaxX,binMinY,binMaxY));
  TH1D *hTemp5 = dynamic_cast<TH1D *>(fHistN[iAnalysisType]->ProjectionZ(Form("%s_Cent_%.0f_%.0f_Psi_%.0f_%.0f",fHistN[iAnalysisType]->GetName(),centrMin,centrMax,psiMin,psiMax),binMinX,binMaxX,binMinY,binMaxY));
  TH1D *hTemp6 = dynamic_cast<TH1D *>(fHistP[iAnalysisType]->ProjectionZ(Form("%s_Cent_%.0f_%.0f_Psi_%.0f_%.0f",fHistP[iAnalysisType]->GetName(),centrMin,centrMax,psiMin,psiMax),binMinX,binMaxX,binMinY,binMaxY));

  if((hTemp1)&&(hTemp2)&&(hTemp3)&&(hTemp4)) {
    hTemp1->Sumw2();
    hTemp2->Sumw2();
    hTemp3->Sumw2();
    hTemp4->Sumw2();
    hTemp1->Add(hTemp3,-2.);
    hTemp1->Scale(1./hTemp5->GetEntries());
    hTemp2->Add(hTemp4,-2.);
    hTemp2->Scale(1./hTemp6->GetEntries());
    gHistBalanceFunctionHistogram->Add(hTemp1,hTemp2,1.,1.);
    gHistBalanceFunctionHistogram->Scale(0.5/fP2Step[iAnalysisType]);
  }

  PrintResults(iAnalysisType,gHistBalanceFunctionHistogram);

  return gHistBalanceFunctionHistogram;
}
