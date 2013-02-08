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
#include <TFile.h>
#include <TF1.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TGraphErrors.h>
#include <TString.h>

#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"

#include "AliBalance.h"

using std::cout;
using std::cerr;
using std::endl;

ClassImp(AliBalance)

//____________________________________________________________________//
AliBalance::AliBalance() :
  TObject(), 
  fShuffle(kFALSE),
  fHBTcut(kFALSE),
  fConversionCut(kFALSE),
  fAnalysisLevel("ESD"),
  fAnalyzedEvents(0) ,
  fCentralityId(0) ,
  fCentStart(0.),
  fCentStop(0.),
  fHistHBTbefore(NULL),
  fHistHBTafter(NULL),
  fHistConversionbefore(NULL),
  fHistConversionafter(NULL)
{
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
  TObject(balance), 
  fShuffle(balance.fShuffle),
  fHBTcut(balance.fHBTcut), 
  fConversionCut(balance.fConversionCut), 
  fAnalysisLevel(balance.fAnalysisLevel),
  fAnalyzedEvents(balance.fAnalyzedEvents), 
  fCentralityId(balance.fCentralityId),
  fCentStart(balance.fCentStart),
  fCentStop(balance.fCentStop),
  fHistHBTbefore(balance.fHistHBTbefore),
  fHistHBTafter(balance.fHistHBTafter),
  fHistConversionbefore(balance.fHistConversionbefore),
  fHistConversionafter(balance.fHistConversionafter) {
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
    fCentStart       = balance.fCentStart;
    fCentStop        = balance.fCentStop; 

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

  // global switch disabling the reference 
  // (to avoid "Replacing existing TH1" if several wagons are created in train)
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  TString histName;
  for(Int_t iAnalysisType = 0; iAnalysisType < ANALYSIS_TYPES; iAnalysisType++) {
    histName = "fHistP"; histName += kBFAnalysisType[iAnalysisType]; 
    if(fShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistP[iAnalysisType] = new TH2D(histName.Data(),"",fCentStop-fCentStart,fCentStart,fCentStop,100,fP1Start[iAnalysisType],fP1Stop[iAnalysisType]);

    histName = "fHistN"; histName += kBFAnalysisType[iAnalysisType]; 
    if(fShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistN[iAnalysisType] = new TH2D(histName.Data(),"",fCentStop-fCentStart,fCentStart,fCentStop,100,fP1Start[iAnalysisType],fP1Stop[iAnalysisType]);
  
    histName = "fHistPN"; histName += kBFAnalysisType[iAnalysisType]; 
    if(fShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistPN[iAnalysisType] = new TH2D(histName.Data(),"",fCentStop-fCentStart,fCentStart,fCentStop,fNumberOfBins[iAnalysisType],fP2Start[iAnalysisType],fP2Stop[iAnalysisType]);
    
    histName = "fHistNP"; histName += kBFAnalysisType[iAnalysisType]; 
    if(fShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistNP[iAnalysisType] = new TH2D(histName.Data(),"",fCentStop-fCentStart,fCentStart,fCentStop,fNumberOfBins[iAnalysisType],fP2Start[iAnalysisType],fP2Stop[iAnalysisType]);
    
    histName = "fHistPP"; histName += kBFAnalysisType[iAnalysisType]; 
    if(fShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistPP[iAnalysisType] = new TH2D(histName.Data(),"",fCentStop-fCentStart,fCentStart,fCentStop,fNumberOfBins[iAnalysisType],fP2Start[iAnalysisType],fP2Stop[iAnalysisType]);
    
    histName = "fHistNN"; histName += kBFAnalysisType[iAnalysisType]; 
    if(fShuffle) histName.Append("_shuffle");
    if(fCentralityId) histName += fCentralityId.Data();
    fHistNN[iAnalysisType] = new TH2D(histName.Data(),"",fCentStop-fCentStart,fCentStart,fCentStop,fNumberOfBins[iAnalysisType],fP2Start[iAnalysisType],fP2Stop[iAnalysisType]);
  }

  // QA histograms
  fHistHBTbefore        = new TH2D("fHistHBTbefore","before HBT cut",200,0,2,200,0,200);
  fHistHBTafter         = new TH2D("fHistHBTafter","after HBT cut",200,0,2,200,0,200);
  fHistConversionbefore = new TH2D("fHistConversionbefore","before Conversion cut",200,0,2,200,0,200);
  fHistConversionafter  = new TH2D("fHistConversionafter","after Conversion cut",200,0,2,200,0,200);

  TH1::AddDirectory(oldStatus);

}

//____________________________________________________________________//
void AliBalance::PrintAnalysisSettings() {
  //prints the analysis settings
  
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
void AliBalance::CalculateBalance(Float_t fCentrality,vector<Double_t> **chargeVector,Float_t bSign) {
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
  //Printf("(AliBalance) Number of tracks: %d",gNtrack);

  for(i = 0; i < gNtrack;i++){
      Short_t charge          = chargeVector[0]->at(i);
      Double_t rapidity       = chargeVector[1]->at(i);
      Double_t pseudorapidity = chargeVector[2]->at(i);
      Double_t phi            = chargeVector[3]->at(i);
      
      //0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi
      for(Int_t iAnalysisType = 0; iAnalysisType < ANALYSIS_TYPES; iAnalysisType++) {
	if(iAnalysisType == kEta) {
	  if((pseudorapidity >= fP1Start[iAnalysisType]) && (pseudorapidity <= fP1Stop[iAnalysisType])) {
	    if(charge > 0) {
	      fNp[iAnalysisType] += 1.;
	      fHistP[iAnalysisType]->Fill(fCentrality,pseudorapidity);
	    }//charge > 0
	    if(charge < 0) {
	      fNn[iAnalysisType] += 1.;
	      fHistN[iAnalysisType]->Fill(fCentrality,pseudorapidity);
	    }//charge < 0
	  }//p1 interval check
	}//analysis type: eta
	else if(iAnalysisType == kPhi) {
	  if((phi >= fP1Start[iAnalysisType]) && (phi <= fP1Stop[iAnalysisType])) {
	    if(charge > 0) {
	      fNp[iAnalysisType] += 1.;
	      fHistP[iAnalysisType]->Fill(fCentrality,phi);
	    }//charge > 0
	    if(charge < 0) {
	      fNn[iAnalysisType] += 1.;
	      fHistN[iAnalysisType]->Fill(fCentrality,phi);
	    }//charge < 0
	  }//p1 interval check
	}//analysis type: phi
	else {
	  if((rapidity >= fP1Start[iAnalysisType]) && (rapidity <= fP1Stop[iAnalysisType])) {
	    if(charge > 0) {
	      fNp[iAnalysisType] += 1.;
	      fHistP[iAnalysisType]->Fill(fCentrality,rapidity);
	    }//charge > 0
	    if(charge < 0) {
	      fNn[iAnalysisType] += 1.;
	      fHistN[iAnalysisType]->Fill(fCentrality,rapidity);
	    }//charge < 0
	  }//p1 interval check
	}//analysis type: y, qside, qout, qlong, qinv
      }//analysis type loop
  }

  //Printf("Np: %lf - Nn: %lf",fNp[0],fNn[0]);

  Double_t dy = 0., deta = 0.;
  Double_t qLong = 0., qOut = 0., qSide = 0., qInv = 0.;
  Double_t dphi = 0.;

  Short_t charge1  = 0;
  Double_t eta1 = 0., rap1 = 0.;
  Double_t px1 = 0., py1 = 0., pz1 = 0.;
  Double_t pt1 = 0.;
  Double_t energy1 = 0.;
  Double_t phi1    = 0.;

  Short_t charge2  = 0;
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
    
    for(j = 0; j < i; j++) {
   
      charge2 = chargeVector[0]->at(j);
      rap2    = chargeVector[1]->at(j);
      eta2    = chargeVector[2]->at(j);
      phi2    = chargeVector[3]->at(j);
      px2     = chargeVector[4]->at(j);
      py2     = chargeVector[5]->at(j);
      pz2     = chargeVector[6]->at(j);
      pt2     = chargeVector[7]->at(j);
      energy2 = chargeVector[8]->at(j);

	// filling the arrays

	// RAPIDITY 
        dy = TMath::Abs(rap1 - rap2);

	// Eta
	deta = TMath::Abs(eta1 - eta2);

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
	dphi = TMath::Abs(phi1 - phi2);
	if(dphi>180) dphi = 360 - dphi;  //dphi should be between 0 and 180!

	// HBT like cut
	if(fHBTcut && charge1 * charge2 > 0){
	  //if( dphi < 3 || deta < 0.01 ){   // VERSION 1
	  //  continue;
	  
	  // VERSION 2 (Taken from DPhiCorrelations)
	  // the variables & cuthave been developed by the HBT group 
	  // see e.g. https://indico.cern.ch/materialDisplay.py?contribId=36&sessionId=6&materialId=slides&confId=142700

	  fHistHBTbefore->Fill(deta,dphi);
	  
	  // optimization
	  if (TMath::Abs(deta) < 0.02 * 2.5 * 3) //twoTrackEfficiencyCutValue = 0.02 [default for dphicorrelations]
	    {

	      // phi in rad
	      Float_t phi1rad = phi1*TMath::DegToRad();
	      Float_t phi2rad = phi2*TMath::DegToRad();

	      // check first boundaries to see if is worth to loop and find the minimum
	      Float_t dphistar1 = GetDPhiStar(phi1rad, pt1, charge1, phi2rad, pt2, charge2, 0.8, bSign);
	      Float_t dphistar2 = GetDPhiStar(phi1rad, pt1, charge1, phi2rad, pt2, charge2, 2.5, bSign);
	      
	      const Float_t kLimit = 0.02 * 3;
	      
	      Float_t dphistarminabs = 1e5;
	      Float_t dphistarmin = 1e5;

	      if (TMath::Abs(dphistar1) < kLimit || TMath::Abs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0 )
		{
		  for (Double_t rad=0.8; rad<2.51; rad+=0.01) 
		    {
		      Float_t dphistar = GetDPhiStar(phi1rad, pt1, charge1, phi2rad, pt2, charge2, rad, bSign);
		      Float_t dphistarabs = TMath::Abs(dphistar);
		      
		      if (dphistarabs < dphistarminabs)
			{
			  dphistarmin = dphistar;
			  dphistarminabs = dphistarabs;
			}
		    }
		
		  if (dphistarminabs < 0.02 && TMath::Abs(deta) < 0.02)
		    {
		      //AliInfo(Form("HBT: Removed track pair %d %d with [[%f %f]] %f %f %f | %f %f %d %f %f %d %f", i, j, deta, dphi, dphistarminabs, dphistar1, dphistar2, phi1rad, pt1, charge1, phi2rad, pt2, charge2, bSign));
		      continue;
		    }
		}
	    }
	  fHistHBTafter->Fill(deta,dphi);
	}
	
	// conversions
	if(fConversionCut){
	  if (charge1 * charge2 < 0)
	    {

	      fHistConversionbefore->Fill(deta,dphi);

	      Float_t m0 = 0.510e-3;
	      Float_t tantheta1 = 1e10;

	      // phi in rad
	      Float_t phi1rad = phi1*TMath::DegToRad();
	      Float_t phi2rad = phi2*TMath::DegToRad();
	      
	      if (eta1 < -1e-10 || eta1 > 1e-10)
		tantheta1 = 2 * TMath::Exp(-eta1) / ( 1 - TMath::Exp(-2*eta1));
	      
	      Float_t tantheta2 = 1e10;
	      if (eta2 < -1e-10 || eta2 > 1e-10)
		tantheta2 = 2 * TMath::Exp(-eta2) / ( 1 - TMath::Exp(-2*eta2));
	      
	      Float_t e1squ = m0 * m0 + pt1 * pt1 * (1.0 + 1.0 / tantheta1 / tantheta1);
	      Float_t e2squ = m0 * m0 + pt2 * pt2 * (1.0 + 1.0 / tantheta2 / tantheta2);

	      Float_t masssqu = 2 * m0 * m0 + 2 * ( TMath::Sqrt(e1squ * e2squ) - ( pt1 * pt2 * ( TMath::Cos(phi1rad - phi2rad) + 1.0 / tantheta1 / tantheta2 ) ) );

	      if (masssqu < 0.04*0.04){
		//AliInfo(Form("Conversion: Removed track pair %d %d with [[%f %f] %f %f] %d %d <- %f %f  %f %f   %f %f ", i, j, deta, dphi, masssqu, charge1, charge2,eta1,eta2,phi1,phi2,pt1,pt2));
		continue;
	      }
	      fHistConversionafter->Fill(deta,dphi);
	    }
	}
	
	
	//0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi
	if((rap1 >= fP1Start[kRapidity]) && (rap1 <= fP1Stop[kRapidity]) && (rap2 >= fP1Start[kRapidity]) && (rap2 <= fP1Stop[kRapidity])) {

	  // rapidity
	  if( dy > fP2Start[kRapidity] && dy < fP2Stop[kRapidity]){
	    iBin = Int_t((dy-fP2Start[kRapidity])/fP2Step[kRapidity]);
	    if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
	      
	      if((charge1 > 0)&&(charge2 > 0)) {
		fNpp[kRapidity][iBin] += 1.;
		fHistPP[kRapidity]->Fill(fCentrality,dy);
	      }
	      else if((charge1 < 0)&&(charge2 < 0)) {
		fNnn[kRapidity][iBin] += 1.;
		fHistNN[kRapidity]->Fill(fCentrality,dy);
	      }
	      else if((charge1 > 0)&&(charge2 < 0)) {
		fNpn[kRapidity][iBin] += 1.;
		fHistPN[kRapidity]->Fill(fCentrality,dy);
	      }
	      else if((charge1 < 0)&&(charge2 > 0)) {
		fNpn[kRapidity][iBin] += 1.;
		    fHistPN[kRapidity]->Fill(fCentrality,dy);
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
		fNpp[kEta][iBin] += 1.;
		fHistPP[kEta]->Fill(fCentrality,deta);
	      }
	      if((charge1 < 0)&&(charge2 < 0)) {
		fNnn[kEta][iBin] += 1.;
	   	    fHistNN[kEta]->Fill(fCentrality,deta);
	      }
	      if((charge1 > 0)&&(charge2 < 0)) {
		fNpn[kEta][iBin] += 1.;
		fHistPN[kEta]->Fill(fCentrality,deta);
	      }
	      if((charge1 < 0)&&(charge2 > 0)) {
		fNpn[kEta][iBin] += 1.;
	   	    fHistPN[kEta]->Fill(fCentrality,deta);
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
		fNpp[kQlong][iBin] += 1.;
		fHistPP[kQlong]->Fill(fCentrality,qLong);
	      }
	      if((charge1 < 0)&&(charge2 < 0)) {
		fNnn[kQlong][iBin] += 1.;
		fHistNN[kQlong]->Fill(fCentrality,qLong);
	      }
	      if((charge1 > 0)&&(charge2 < 0)) {
		fNpn[kQlong][iBin] += 1.;
		fHistPN[kQlong]->Fill(fCentrality,qLong);
	      }
	      if((charge1 < 0)&&(charge2 > 0)) {
		fNpn[kQlong][iBin] += 1.;
		fHistPN[kQlong]->Fill(fCentrality,qLong);
	      }
	    }//BF binning check
	  }//p2 interval check
	}//p1 interval check
	  
	if((rap1 >= fP1Start[kRapidity]) && (rap1 <= fP1Stop[kRapidity]) && (rap2 >= fP1Start[kRapidity]) && (rap2 <= fP1Stop[kRapidity])) {
	  if( qOut > fP2Start[kQout] && qOut < fP2Stop[kQout]){
	    iBin = Int_t((qOut-fP2Start[kQout])/fP2Step[kQout]);	
	    if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
	      if((charge1 > 0)&&(charge2 > 0)) {
		fNpp[kQout][iBin] += 1.;
		fHistPP[kQout]->Fill(fCentrality,qOut);
	  	  }
	      if((charge1 < 0)&&(charge2 < 0)) {
		fNnn[kQout][iBin] += 1.;
		fHistNN[kQout]->Fill(fCentrality,qOut);
	      }
	      if((charge1 > 0)&&(charge2 < 0)) {
		fNpn[kQout][iBin] += 1.;
		fHistPN[kQout]->Fill(fCentrality,qOut);
	      }
	      if((charge1 < 0)&&(charge2 > 0)) {
		fNpn[kQout][iBin] += 1.;
		fHistPN[kQout]->Fill(fCentrality,qOut);
	      }
	    }//BF binning check
	  }//p2 interval check
	}//p1 interval check	
	
	if((rap1 >= fP1Start[kRapidity]) && (rap1 <= fP1Stop[kRapidity]) && (rap2 >= fP1Start[kRapidity]) && (rap2 <= fP1Stop[kRapidity])) {
	  if( qSide > fP2Start[kQside] && qSide < fP2Stop[kQside]){
	    iBin = Int_t((qSide-fP2Start[kQside])/fP2Step[kQside]);	
	    if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
	      if((charge1 > 0)&&(charge2 > 0)) {
		fNpp[kQside][iBin] += 1.;
		fHistPP[kQside]->Fill(fCentrality,qSide);
	      }
	      if((charge1 < 0)&&(charge2 < 0)) {
		fNnn[kQside][iBin] += 1.;
		fHistNN[kQside]->Fill(fCentrality,qSide);
	      }
	      if((charge1 > 0)&&(charge2 < 0)) {
		fNpn[kQside][iBin] += 1.;
		fHistPN[kQside]->Fill(fCentrality,qSide);
	      }
	      if((charge1 < 0)&&(charge2 > 0)) {
		fNpn[kQside][iBin] += 1.;
		fHistPN[kQside]->Fill(fCentrality,qSide);
		  	  }
	    }//BF binning check
	  }//p2 interval check
	}//p1 interval check
	
	if((rap1 >= fP1Start[kRapidity]) && (rap1 <= fP1Stop[kRapidity]) && (rap2 >= fP1Start[kRapidity]) && (rap2 <= fP1Stop[kRapidity])) {
	  if( qInv > fP2Start[kQinv] && qInv < fP2Stop[kQinv]){
	    iBin = Int_t((qInv-fP2Start[kQinv])/fP2Step[kQinv]);	
	    if(iBin >=0 && iBin < MAXIMUM_NUMBER_OF_STEPS){
	      if((charge1 > 0)&&(charge2 > 0)) {
		fNpp[kQinv][iBin] += 1.;
		fHistPP[kQinv]->Fill(fCentrality,qInv);
	      }
	      if((charge1 < 0)&&(charge2 < 0)) {
		fNnn[kQinv][iBin] += 1.;
		fHistNN[kQinv]->Fill(fCentrality,qInv);
	      }
	      if((charge1 > 0)&&(charge2 < 0)) {
		fNpn[kQinv][iBin] += 1.;
		fHistPN[kQinv]->Fill(fCentrality,qInv);
	      }
	      if((charge1 < 0)&&(charge2 > 0)) {
		fNpn[kQinv][iBin] += 1.;
		fHistPN[kQinv]->Fill(fCentrality,qInv);
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
		fNpp[kPhi][iBin] += 1.;
		fHistPP[kPhi]->Fill(fCentrality,dphi);
	      }
	      if((charge1 < 0)&&(charge2 < 0)) {
		fNnn[kPhi][iBin] += 1.;
		fHistNN[kPhi]->Fill(fCentrality,dphi);
	      }
	      if((charge1 > 0)&&(charge2 < 0)) {
		fNpn[kPhi][iBin] += 1.;
		fHistPN[kPhi]->Fill(fCentrality,dphi);
	      }
	      if((charge1 < 0)&&(charge2 > 0)) {
		fNpn[kPhi][iBin] += 1.;
		fHistPN[kPhi]->Fill(fCentrality,dphi);
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
  Double_t gSumXi = 0.0, gSumBi = 0.0, gSumBiXi = 0.0;
  Double_t gSumBiXi2 = 0.0, gSumBi2Xi2 = 0.0;
  Double_t gSumDeltaBi2 = 0.0, gSumXi2DeltaBi2 = 0.0;
  Double_t deltaBalP2 = 0.0, integral = 0.0;
  Double_t deltaErrorNew = 0.0;
  
  // cout<<"=================================================="<<endl;
  // for(Int_t i = 1; i <= fNumberOfBins[iAnalysisType]; i++) { 
  //   x[i-1] = fP2Start[iAnalysisType] + fP2Step[iAnalysisType]*i + fP2Step[iAnalysisType]/2;
  //   cout<<"B: "<<gHistBalance->GetBinContent(i)<<"\t Error: "<<gHistBalance->GetBinError(i)<<"\t bin: "<<gHistBalance->GetBinCenter(i)<<endl;
  // } 
  // cout<<"=================================================="<<endl;
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
  integralError *= 1.0;

  Double_t delta = gSumBiXi / gSumBi; delta *= 1.0;
  Double_t deltaError = (gSumBiXi / gSumBi) * TMath::Sqrt(TMath::Power((TMath::Sqrt(gSumXi2DeltaBi2)/gSumBiXi),2) + TMath::Power((gSumDeltaBi2/gSumBi),2) );
  deltaError *= 1.0;
  // cout<<"Analysis type: "<<kBFAnalysisType[iAnalysisType].Data()<<endl;
  // cout<<"Width: "<<delta<<"\t Error: "<<deltaError<<endl;
  // cout<<"New error: "<<deltaErrorNew<<endl;
  // cout<<"Integral: "<<integral<<"\t Error: "<<integralError<<endl;
  // cout<<"=================================================="<<endl;
}
 
//____________________________________________________________________//
TH1D *AliBalance::GetBalanceFunctionHistogram(Int_t iAnalysisType,Double_t centrMin, Double_t centrMax, Double_t etaWindow,Bool_t correctWithEfficiency, Bool_t correctWithAcceptanceOnly, Bool_t correctWithMixed, TH1D *hMixed[4]) {
  //Returns the BF histogram, extracted from the 6 TH2D objects 
  //(private members) of the AliBalance class.
  //
  // Acceptance correction: 
  // - only for analysis type = kEta
  // - only if etaWindow > 0 (default = -1.)
  // - calculated as proposed by STAR 
  //
  TString gAnalysisType[ANALYSIS_TYPES] = {"y","eta","qlong","qout","qside","qinv","phi"};
  TString histName = "gHistBalanceFunctionHistogram";
  histName += gAnalysisType[iAnalysisType];

  SetInterval(iAnalysisType, fHistP[iAnalysisType]->GetYaxis()->GetXmin(),
	      fHistP[iAnalysisType]->GetYaxis()->GetXmin(),
	      fHistPP[iAnalysisType]->GetNbinsY(),
	      fHistPP[iAnalysisType]->GetYaxis()->GetXmin(),
	      fHistPP[iAnalysisType]->GetYaxis()->GetXmax());

  // determine the projection thresholds
  Int_t binMinX, binMinY, binMinZ;
  Int_t binMaxX, binMaxY, binMaxZ;

  fHistPP[iAnalysisType]->GetBinXYZ(fHistPP[iAnalysisType]->FindBin(centrMin),binMinX,binMinY,binMinZ);
  fHistPP[iAnalysisType]->GetBinXYZ(fHistPP[iAnalysisType]->FindBin(centrMax),binMaxX,binMaxY,binMaxZ);

  TH1D *gHistBalanceFunctionHistogram = new TH1D(histName.Data(),"",fHistPP[iAnalysisType]->GetNbinsY(),fHistPP[iAnalysisType]->GetYaxis()->GetXmin(),fHistPP[iAnalysisType]->GetYaxis()->GetXmax());
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

  TH1D *hTemp1 = dynamic_cast<TH1D *>(fHistPN[iAnalysisType]->ProjectionY(Form("%s_Cent_%.0f_%.0f",fHistPN[iAnalysisType]->GetName(),centrMin,centrMax),binMinX,binMaxX));
  TH1D *hTemp2 = dynamic_cast<TH1D *>(fHistPN[iAnalysisType]->ProjectionY(Form("%s_Cent_%.0f_%.0f_copy",fHistPN[iAnalysisType]->GetName(),centrMin,centrMax),binMinX,binMaxX));
  TH1D *hTemp3 = dynamic_cast<TH1D *>(fHistNN[iAnalysisType]->ProjectionY(Form("%s_Cent_%.0f_%.0f",fHistNN[iAnalysisType]->GetName(),centrMin,centrMax),binMinX,binMaxX));
  TH1D *hTemp4 = dynamic_cast<TH1D *>(fHistPP[iAnalysisType]->ProjectionY(Form("%s_Cent_%.0f_%.0f",fHistPP[iAnalysisType]->GetName(),centrMin,centrMax),binMinX,binMaxX));
  TH1D *hTemp5 = dynamic_cast<TH1D *>(fHistN[iAnalysisType]->ProjectionY(Form("%s_Cent_%.0f_%.0f",fHistN[iAnalysisType]->GetName(),centrMin,centrMax),binMinX,binMaxX));
  TH1D *hTemp6 = dynamic_cast<TH1D *>(fHistP[iAnalysisType]->ProjectionY(Form("%s_Cent_%.0f_%.0f",fHistP[iAnalysisType]->GetName(),centrMin,centrMax),binMinX,binMaxX));

  // get the file with the efficiency matrices
  // withAcceptanceOnly: Data single distributions are normalized to 1 (efficiency not taken into account)
  // else : Data single distributions are normalized to give single particle efficiency of MC
  TFile *fEfficiencyMatrix = NULL;
  if(correctWithEfficiency || correctWithMixed){
    if(correctWithAcceptanceOnly) fEfficiencyMatrix = TFile::Open("$ALICE_ROOT/PWGCF/EBYE/macros/accOnlyFromConvolutionAllCent.root");
    else  fEfficiencyMatrix = TFile::Open("$ALICE_ROOT/PWGCF/EBYE/macros/effFromConvolutionAllCent.root");
    if(!fEfficiencyMatrix){
      AliError("Efficiency histogram file not found");
      return NULL;
    }
  }

  // do correction with the efficiency calculated from MC + Data (for single particles and two particle correlations)
  // - single particle efficiencies from MC (AliAnalysiTaskEfficiency)
  // - two particle efficiencies from convolution of data single particle distributions 
  //   (normalized to single particle efficiency)
  if(iAnalysisType == kEta && etaWindow > 0 && correctWithEfficiency && !correctWithMixed){

    TH1F* hEffP  = NULL;
    TH1F* hEffN  = NULL;
    TH1F* hEffPP = (TH1F*)fEfficiencyMatrix->Get(Form("etaEffPP_Cent%.0f-%.0f_Data",centrMin,centrMax));
    TH1F* hEffNN = (TH1F*)fEfficiencyMatrix->Get(Form("etaEffNN_Cent%.0f-%.0f_Data",centrMin,centrMax));
    TH1F* hEffPN = (TH1F*)fEfficiencyMatrix->Get(Form("etaEffPN_Cent%.0f-%.0f_Data",centrMin,centrMax));

    // take the data distributions
    if(correctWithAcceptanceOnly){
      hEffP = (TH1F*)fEfficiencyMatrix->Get(Form("etaEffP_Cent%.0f-%.0f_Data",centrMin,centrMax));
      hEffN = (TH1F*)fEfficiencyMatrix->Get(Form("etaEffN_Cent%.0f-%.0f_Data",centrMin,centrMax));
    }
    // take the MC distributions
    else{
      hEffP = (TH1F*)fEfficiencyMatrix->Get(Form("etaEffP_Cent%.0f-%.0f_MC",centrMin,centrMax));
      hEffN = (TH1F*)fEfficiencyMatrix->Get(Form("etaEffN_Cent%.0f-%.0f_MC",centrMin,centrMax));
    }

    if( !hEffP || !hEffN || !hEffPP || !hEffNN || !hEffPN){
      AliError(Form("Efficiency (eta) histograms not found: etaEffPP_Cent%.0f-%.0f_Data",centrMin,centrMax));
      return NULL;
    }

    for(Int_t iBin = 0; iBin < hEffP->GetNbinsX(); iBin++){
      hTemp5->SetBinError(iBin+1,hTemp5->GetBinError(iBin+1)/hEffN->GetBinContent(hEffN->FindBin(hTemp5->GetBinCenter(iBin+1))));
      hTemp5->SetBinContent(iBin+1,hTemp5->GetBinContent(iBin+1)/hEffN->GetBinContent(hEffN->FindBin(hTemp5->GetBinCenter(iBin+1))));

      hTemp6->SetBinError(iBin+1,hTemp6->GetBinError(iBin+1)/hEffP->GetBinContent(hEffP->FindBin(hTemp6->GetBinCenter(iBin+1))));
      hTemp6->SetBinContent(iBin+1,hTemp6->GetBinContent(iBin+1)/hEffP->GetBinContent(hEffP->FindBin(hTemp6->GetBinCenter(iBin+1))));
    }
  
    for(Int_t iBin = 0; iBin < gHistBalanceFunctionHistogram->GetNbinsX(); iBin++){

      hTemp1->SetBinError(iBin+1,hTemp1->GetBinError(iBin+1)/hEffPN->GetBinContent(hEffPN->FindBin(hTemp1->GetBinCenter(iBin+1))));
      hTemp1->SetBinContent(iBin+1,hTemp1->GetBinContent(iBin+1)/hEffPN->GetBinContent(hEffPN->FindBin(hTemp1->GetBinCenter(iBin+1))));
      hTemp2->SetBinError(iBin+1,hTemp2->GetBinError(iBin+1)/hEffPN->GetBinContent(hEffPN->FindBin(hTemp2->GetBinCenter(iBin+1))));
      hTemp2->SetBinContent(iBin+1,hTemp2->GetBinContent(iBin+1)/hEffPN->GetBinContent(hEffPN->FindBin(hTemp2->GetBinCenter(iBin+1))));
      hTemp3->SetBinError(iBin+1,hTemp3->GetBinError(iBin+1)/hEffNN->GetBinContent(hEffNN->FindBin(hTemp3->GetBinCenter(iBin+1))));
      hTemp3->SetBinContent(iBin+1,hTemp3->GetBinContent(iBin+1)/hEffNN->GetBinContent(hEffNN->FindBin(hTemp3->GetBinCenter(iBin+1))));
      hTemp4->SetBinError(iBin+1,hTemp4->GetBinError(iBin+1)/hEffPP->GetBinContent(hEffPP->FindBin(hTemp4->GetBinCenter(iBin+1))));      
      hTemp4->SetBinContent(iBin+1,hTemp4->GetBinContent(iBin+1)/hEffPP->GetBinContent(hEffPP->FindBin(hTemp4->GetBinCenter(iBin+1))));
      
    }

    // TF1 *fPP = new TF1("fPP","pol1",0,1.6);  // phase space factor + efficiency for ++
    // fPP->SetParameters(0.736466,-0.461529);
    // TF1 *fNN = new TF1("fNN","pol1",0,1.6);  // phase space factor + efficiency for --
    // fNN->SetParameters(0.718616,-0.450473);
    // TF1 *fPN = new TF1("fPN","pol1",0,1.6);  // phase space factor + efficiency for +-
    // fPN->SetParameters(0.727507,-0.455981);
    
    // for(Int_t iBin = 0; iBin < gHistBalanceFunctionHistogram->GetNbinsX(); iBin++){
    //   hTemp1->SetBinContent(iBin+1,hTemp1->GetBinContent(iBin+1)/fPN->Eval(hTemp1->GetBinCenter(iBin+1)));
    //   hTemp1->SetBinError(iBin+1,hTemp1->GetBinError(iBin+1)/fPN->Eval(hTemp1->GetBinCenter(iBin+1)));
    //   hTemp2->SetBinContent(iBin+1,hTemp2->GetBinContent(iBin+1)/fPN->Eval(hTemp1->GetBinCenter(iBin+1)));
    //   hTemp2->SetBinError(iBin+1,hTemp2->GetBinError(iBin+1)/fPN->Eval(hTemp1->GetBinCenter(iBin+1)));
    //   hTemp3->SetBinContent(iBin+1,hTemp3->GetBinContent(iBin+1)/fNN->Eval(hTemp1->GetBinCenter(iBin+1)));
    //   hTemp3->SetBinError(iBin+1,hTemp3->GetBinError(iBin+1)/fNN->Eval(hTemp1->GetBinCenter(iBin+1)));
    //   hTemp4->SetBinContent(iBin+1,hTemp4->GetBinContent(iBin+1)/fPP->Eval(hTemp1->GetBinCenter(iBin+1)));
    //   hTemp4->SetBinError(iBin+1,hTemp4->GetBinError(iBin+1)/fPP->Eval(hTemp1->GetBinCenter(iBin+1)));
    // }      
  }

  // do correction with the efficiency calculated from MC + Data (for single particles and two particle correlations)
  // - single particle efficiencies from MC (AliAnalysiTaskEfficiency)
  // - two particle efficiencies from convolution of data single particle distributions 
  //   (normalized to single particle efficiency)  
  if(iAnalysisType == kPhi && correctWithEfficiency && !correctWithMixed){

    TH1F* hEffPhiP  = NULL;
    TH1F* hEffPhiN  = NULL;
    TH1F* hEffPhiPP = (TH1F*)fEfficiencyMatrix->Get(Form("phiEffPP_Cent%.0f-%.0f_Data",centrMin,centrMax));
    TH1F* hEffPhiNN = (TH1F*)fEfficiencyMatrix->Get(Form("phiEffNN_Cent%.0f-%.0f_Data",centrMin,centrMax));
    TH1F* hEffPhiPN = (TH1F*)fEfficiencyMatrix->Get(Form("phiEffPN_Cent%.0f-%.0f_Data",centrMin,centrMax));

    // take the data distributions
    if(correctWithAcceptanceOnly){
      hEffPhiP = (TH1F*)fEfficiencyMatrix->Get(Form("phiEffP_Cent%.0f-%.0f_Data",centrMin,centrMax));
      hEffPhiN = (TH1F*)fEfficiencyMatrix->Get(Form("phiEffN_Cent%.0f-%.0f_Data",centrMin,centrMax));
    }
    // take the MC distributions
    else{
      hEffPhiP = (TH1F*)fEfficiencyMatrix->Get(Form("phiEffP_Cent%.0f-%.0f_MC",centrMin,centrMax));
      hEffPhiN = (TH1F*)fEfficiencyMatrix->Get(Form("phiEffN_Cent%.0f-%.0f_MC",centrMin,centrMax));
    }

    if( !hEffPhiP || !hEffPhiN || !hEffPhiPP || !hEffPhiNN || !hEffPhiPN){
      AliError("Efficiency (phi) histograms not found");
      return NULL;
    }

    for(Int_t iBin = 0; iBin < hEffPhiP->GetNbinsX(); iBin++){
      hTemp5->SetBinError(iBin+1,hTemp5->GetBinError(iBin+1)/hEffPhiN->GetBinContent(hEffPhiN->FindBin(hTemp5->GetBinCenter(iBin+1))));
      hTemp5->SetBinContent(iBin+1,hTemp5->GetBinContent(iBin+1)/hEffPhiN->GetBinContent(hEffPhiN->FindBin(hTemp5->GetBinCenter(iBin+1))));

      hTemp6->SetBinError(iBin+1,hTemp6->GetBinError(iBin+1)/hEffPhiP->GetBinContent(hEffPhiP->FindBin(hTemp6->GetBinCenter(iBin+1))));
      hTemp6->SetBinContent(iBin+1,hTemp6->GetBinContent(iBin+1)/hEffPhiP->GetBinContent(hEffPhiP->FindBin(hTemp6->GetBinCenter(iBin+1))));
    }

    for(Int_t iBin = 0; iBin < gHistBalanceFunctionHistogram->GetNbinsX(); iBin++){

      hTemp1->SetBinError(iBin+1,hTemp1->GetBinError(iBin+1)/hEffPhiPN->GetBinContent(hEffPhiPN->FindBin(hTemp1->GetBinCenter(iBin+1))));
      hTemp1->SetBinContent(iBin+1,hTemp1->GetBinContent(iBin+1)/hEffPhiPN->GetBinContent(hEffPhiPN->FindBin(hTemp1->GetBinCenter(iBin+1))));
      hTemp2->SetBinError(iBin+1,hTemp2->GetBinError(iBin+1)/hEffPhiPN->GetBinContent(hEffPhiPN->FindBin(hTemp2->GetBinCenter(iBin+1))));
      hTemp2->SetBinContent(iBin+1,hTemp2->GetBinContent(iBin+1)/hEffPhiPN->GetBinContent(hEffPhiPN->FindBin(hTemp2->GetBinCenter(iBin+1))));
      hTemp3->SetBinError(iBin+1,hTemp3->GetBinError(iBin+1)/hEffPhiNN->GetBinContent(hEffPhiNN->FindBin(hTemp3->GetBinCenter(iBin+1))));
      hTemp3->SetBinContent(iBin+1,hTemp3->GetBinContent(iBin+1)/hEffPhiNN->GetBinContent(hEffPhiNN->FindBin(hTemp3->GetBinCenter(iBin+1))));
      hTemp4->SetBinError(iBin+1,hTemp4->GetBinError(iBin+1)/hEffPhiPP->GetBinContent(hEffPhiPP->FindBin(hTemp4->GetBinCenter(iBin+1))));      
      hTemp4->SetBinContent(iBin+1,hTemp4->GetBinContent(iBin+1)/hEffPhiPP->GetBinContent(hEffPhiPP->FindBin(hTemp4->GetBinCenter(iBin+1))));
      
    }  
  }

  // do the correction with the event mixing directly!
  if(correctWithMixed){
  
    // take the MC distributions (for average efficiency)
    TH1F* hEffP = (TH1F*)fEfficiencyMatrix->Get(Form("etaEffP_Cent%.0f-%.0f_MC",centrMin,centrMax));
    TH1F* hEffN = (TH1F*)fEfficiencyMatrix->Get(Form("etaEffN_Cent%.0f-%.0f_MC",centrMin,centrMax));  

    TH1F* hEffPP = (TH1F*)fEfficiencyMatrix->Get(Form("phiEffPP_Cent%.0f-%.0f_Data",centrMin,centrMax));
    TH1F* hEffNN = (TH1F*)fEfficiencyMatrix->Get(Form("phiEffNN_Cent%.0f-%.0f_Data",centrMin,centrMax));
    TH1F* hEffPN = (TH1F*)fEfficiencyMatrix->Get(Form("phiEffPN_Cent%.0f-%.0f_Data",centrMin,centrMax));

    if( !hEffP || !hEffN){
      AliError(Form("Efficiency (eta) histograms not found: etaEffPP_Cent%.0f-%.0f_Data",centrMin,centrMax));
      return NULL;
    }

    if(hMixed[0] && hMixed[1] && hMixed[2] && hMixed[3]){
    
      // scale to average efficiency in the pt region (0.3-1.5) and |eta| < 0.8
      // by multiplying the average single particle efficiencies from HIJING
      // here we assume that the distributions are 1:
      // - in the integral for dphi (for averaging over sector structure)
      // - in the maximum for deta
      Double_t normPMC = (Double_t)hEffP->Integral()/(Double_t)hEffP->GetNbinsX();
      Double_t normNMC = (Double_t)hEffN->Integral()/(Double_t)hEffN->GetNbinsX();
      Double_t normPPMC = (Double_t)hEffPP->Integral()/(Double_t)hEffPP->GetNbinsX();
      Double_t normNNMC = (Double_t)hEffNN->Integral()/(Double_t)hEffNN->GetNbinsX();
      Double_t normPNMC = (Double_t)hEffPN->Integral()/(Double_t)hEffPN->GetNbinsX();

      hMixed[0]->Scale(normPNMC);
      hMixed[1]->Scale(normPNMC);
      hMixed[2]->Scale(normNNMC);
      hMixed[3]->Scale(normPPMC);

      // divide by event mixing
      hTemp1->Divide(hMixed[0]);
      hTemp2->Divide(hMixed[1]);
      hTemp3->Divide(hMixed[2]);
      hTemp4->Divide(hMixed[3]);

      // scale also single histograms with average efficiency
      hTemp5->Scale(1./normNMC);
      hTemp6->Scale(1./normPMC);

    }
    else{
      AliError("Correction with EventMixing requested, but not all Histograms there!");
      return NULL;
    }
  }


  if((hTemp1)&&(hTemp2)&&(hTemp3)&&(hTemp4)) {
    hTemp1->Sumw2();
    hTemp2->Sumw2();
    hTemp3->Sumw2();
    hTemp4->Sumw2();
    hTemp1->Add(hTemp3,-2.);
    hTemp1->Scale(1./hTemp5->Integral());
    hTemp2->Add(hTemp4,-2.);
    hTemp2->Scale(1./hTemp6->Integral());
    gHistBalanceFunctionHistogram->Add(hTemp1,hTemp2,1.,1.);
    gHistBalanceFunctionHistogram->Scale(0.5/fP2Step[iAnalysisType]);
  }

  // do the acceptance correction (only for Eta and etaWindow > 0)
  if(iAnalysisType == kEta && etaWindow > 0 && !correctWithEfficiency && !correctWithMixed){
    for(Int_t iBin = 0; iBin < gHistBalanceFunctionHistogram->GetNbinsX(); iBin++){
      
      Double_t notCorrected = gHistBalanceFunctionHistogram->GetBinContent(iBin+1);
      Double_t corrected    = notCorrected / (1 - (gHistBalanceFunctionHistogram->GetBinCenter(iBin+1))/ etaWindow );
      gHistBalanceFunctionHistogram->SetBinContent(iBin+1, corrected);
      gHistBalanceFunctionHistogram->SetBinError(iBin+1,corrected/notCorrected*gHistBalanceFunctionHistogram->GetBinError(iBin+1));
      
    }
  }
  
  if(fEfficiencyMatrix)   fEfficiencyMatrix->Close();

  PrintResults(iAnalysisType,gHistBalanceFunctionHistogram);

  return gHistBalanceFunctionHistogram;
}
