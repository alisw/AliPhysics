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
  fAnalysisLevel("ESD"),
  fAnalyzedEvents(0) {
  // Default constructor
 
  for(Int_t i = 0; i < ANALYSIS_TYPES; i++){
    if(i == 6) {
      fNumberOfBins[i] = 180;
      fP1Start[i]      = 0.0;
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
      fB[i][j] = 0.0;
      ferror[i][j] = 0.0;
    } 
  }

  TString gAnalysisType[ANALYSIS_TYPES] = {"y","eta","qlong","qout","qside","qinv","phi"};
  TString histName;
  for(Int_t iAnalysisType = 0; iAnalysisType < ANALYSIS_TYPES; iAnalysisType++) {
    histName = "fHistP"; histName += gAnalysisType[iAnalysisType];
    fHistP[iAnalysisType] = new TH1D(histName.Data(),"",100,fP1Start[iAnalysisType],fP1Stop[iAnalysisType]);
    histName = "fHistN"; histName += gAnalysisType[iAnalysisType];
    fHistN[iAnalysisType] = new TH1D(histName.Data(),"",100,fP1Start[iAnalysisType],fP1Stop[iAnalysisType]);
  
    histName = "fHistPN"; histName += gAnalysisType[iAnalysisType];
    fHistPN[iAnalysisType] = new TH1D(histName.Data(),"",fNumberOfBins[iAnalysisType],fP2Start[iAnalysisType],fP2Stop[iAnalysisType]);
    histName = "fHistNP"; histName += gAnalysisType[iAnalysisType];
    fHistNP[iAnalysisType] = new TH1D(histName.Data(),"",fNumberOfBins[iAnalysisType],fP2Start[iAnalysisType],fP2Stop[iAnalysisType]);
    histName = "fHistPP"; histName += gAnalysisType[iAnalysisType];
    fHistPP[iAnalysisType] = new TH1D(histName.Data(),"",fNumberOfBins[iAnalysisType],fP2Start[iAnalysisType],fP2Stop[iAnalysisType]);
    histName = "fHistNN"; histName += gAnalysisType[iAnalysisType];
    fHistNN[iAnalysisType] = new TH1D(histName.Data(),"",fNumberOfBins[iAnalysisType],fP2Start[iAnalysisType],fP2Stop[iAnalysisType]);

  }
}


//____________________________________________________________________//
AliBalance::AliBalance(const AliBalance& balance):
  TObject(balance), fAnalysisLevel(balance.fAnalysisLevel),
  fAnalyzedEvents(balance.fAnalyzedEvents) {
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
      fB[i][j] = 0.0;
      ferror[i][j] = 0.0;
    } 
  }
 }
 

//____________________________________________________________________//
AliBalance::~AliBalance() {
  // Destructor

}


//____________________________________________________________________//
void AliBalance::SetNumberOfBins(Int_t ibin, Int_t ibins) {
  // Sets the number of bins for the analyzed interval
  // Set the same Information for all analyses
  if(ibin == -1){             
    for(Int_t i = 0; i < ANALYSIS_TYPES; i++){
      fNumberOfBins[i] = ibins;
    }
  }
  // Set the Information for one analysis
  else if(ibin > -1 && ibin < ANALYSIS_TYPES){
    fNumberOfBins[ibin] = ibins;
  }
  else{
    AliError("Wrong ANALYSIS number!");
  }
}

//____________________________________________________________________//
void AliBalance::SetInterval(Double_t p1Start, Double_t p1Stop,
			     Int_t ibin, Double_t p2Start, Double_t p2Stop) {
  // Sets the analyzed interval. 
  // Set the same Information for all analyses
  if(ibin == -1){             
    for(Int_t i = 0; i < ANALYSIS_TYPES; i++){
      fP1Start[i] = p1Start;
      fP1Stop[i] = p1Stop;
      fP2Start[i] = p2Start;
      fP2Stop[i] = p2Stop;
      fP2Step[i] = TMath::Abs(p2Start - p2Stop) / (Double_t)fNumberOfBins[i];
    }
  }
  // Set the Information for one analysis
  else if(ibin > -1 && ibin < ANALYSIS_TYPES){
      fP1Start[ibin] = p1Start;
      fP1Stop[ibin] = p1Stop;
      fP2Start[ibin] = p2Start;
      fP2Stop[ibin] = p2Stop;
      fP2Step[ibin] = TMath::Abs(p2Start - p2Stop) / (Double_t)fNumberOfBins[ibin];
  }
  else{
    AliError("Wrong ANALYSIS number!");
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
void AliBalance::CalculateBalance(TObjArray *gTrackArray) {
  // Calculates the balance function
  fAnalyzedEvents++;
  Int_t i = 0 , j = 0;
  Int_t iBin = 0;
  
  AliVParticle* track = 0;
  AliVParticle* track1 = 0;
  AliVParticle* track2 = 0;
    
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
      Short_t charge  = track->Charge();
      Double_t pseudorapidity = track->Eta();
      Double_t rapidity = track->Y();
      Double_t phi = TMath::ATan(track->Py()/track->Px())*180.0/TMath::Pi();
      
      //0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi
      for(Int_t iAnalysisType = 0; iAnalysisType < ANALYSIS_TYPES; iAnalysisType++) {
	if(iAnalysisType == 1) {
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
	if(iAnalysisType == 6) {
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
  //Printf("Np: %lf - Nn: %lf",fNp,fNn);

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
      charge1 = track1->Charge();
      p1      = track1->P();
      pX1     = track1->Px();
      pY1     = track1->Py();
      pZ1     = track1->Pz();
      energy1 = TMath::Sqrt(TMath::Power(track1->P(),2) +
			    TMath::Power(track1->M(),2));
      phi1    = track1->Phi();

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
	charge2 = track2->Charge();
	p2      = track2->P();
	pX2     = track2->Px();
	pY2     = track2->Py();
	pZ2     = track2->Pz();
	energy2 = TMath::Sqrt(TMath::Power(track2->P(),2) +
			      TMath::Power(track2->M(),2));
	phi2    = track2->Phi();

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
	dphi = phi1 - phi2;

	//0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi
	for(Int_t iAnalysisType = 0; iAnalysisType < ANALYSIS_TYPES; iAnalysisType++) {
	  if(iAnalysisType == kRapidity) {
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
	      }
	    }
	  }//rapidity
	
	  if(iAnalysisType == kEta) {
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
	      }
	    }
	  }//pseudorapidity
	  
	  // Qlong, out, side, inv
	  // thresholds missing!
	  if(iAnalysisType == kQlong) {
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
	      }
	    }
	  }//qLong

	  if(iAnalysisType == kQout) {
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
	      }
	    }
	  }//qOut

	  if(iAnalysisType == kQside) {
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
	      }
	    }
	  }//qSide
	
	  if(iAnalysisType == kQinv) {
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
	      }
	    }
	  }//qInv

	  // Phi
	  if(iAnalysisType == kPhi) {
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
	      }
	    }
	  }//phi
	}//analysis type loop
      }//track2 valid
      else continue;
    }//end of 2nd particle loop
  }//end of 1st particle loop
  //Printf("Number of analyzed events: %i",fAnalyzedEvents);
}  



//____________________________________________________________________//
Double_t AliBalance::GetBalance(Int_t a, Int_t p2) {
  // Returns the value of the balance function in bin p2
  fB[a][p2] = 0.5*(((fNpn[a][p2] - 2.0*fNnn[a][p2])/fNn[a]) + ((fNpn[a][p2] - 2.0*fNpp[a][p2])/fNp[a]))/fP2Step[a];
  
  return fB[a][p2];
}
    
//____________________________________________________________________//
Double_t AliBalance::GetError(Int_t a, Int_t p2) {
  // Returns the error on the BF value for bin p2
  ferror[a][p2] = TMath::Sqrt( Double_t(fNpp[a][p2])/(Double_t(fNp[a])*Double_t(fNp[a])) + Double_t(fNnn[a][p2])/(Double_t(fNn[a])*Double_t(fNn[a])) + Double_t(fNpn[a][p2])*TMath::Power((0.5/Double_t(fNp[a]) + 0.5/Double_t(fNn[a])),2))/fP2Step[a];

  return ferror[a][p2];
}


