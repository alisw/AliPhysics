/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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

/*
   Class AliSpherocityUtils  
   It allows to calculate Spherocity using the official track cuts

 */
//           Please report any Bugs to:
//*************************************************************************
//   author: Antonio Ortiz Velasquez  ( antonio.ortiz@nucleares.unam.mx )
/**************************************************************************


  First version, August 29 2018

 ***************************************************************************/


#include "AliStack.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliVTrack.h"
#include "AliVEvent.h"
#include <TMatrixDSym.h>
#include <TMath.h>
#include <TParticlePDG.h>
#include <TParticle.h>
#include "AliESDUtils.h"
#include "AliESDtrackCuts.h"
#include "AliSpherocityUtils.h"
#include <TFile.h>
#include "AliAODHeader.h"
// STL includes
#include <iostream>
using namespace std;


ClassImp(AliSpherocityUtils)

	//______________________________________________________________________
AliSpherocityUtils::AliSpherocityUtils():TObject(),
	fESDEvent(0),
	fAODEvent(0),
	fMCStack(0),
	fNrec(0),
	fAnalysisMC(kFALSE),
	fTrackFilterGlobal(0),
	fAODFilterGlobal(0),
	fMinMult(3),
	fSizeStep(0.1),
	fIsAbsEta(kTRUE),
	fEtaMaxCut(0.8),
	fEtaMinCut(0.0),
	fPtMaxCut(1e8),
	fPtMinCut(0.15),
	fRunNumber(0)
	
{
	// Default contructor
}


//_____________________________________________________________________________
AliSpherocityUtils::AliSpherocityUtils(const AliSpherocityUtils &c) : 
	TObject(c),
	fESDEvent(0),
	fAODEvent(0),
	fMCStack(0),
	fNrec(0),
	fAnalysisMC(kFALSE),
	fTrackFilterGlobal(0),
	fAODFilterGlobal(0),
	fMinMult(3),
	fSizeStep(0.1),
	fIsAbsEta(kTRUE),
	fEtaMaxCut(0.8),
	fEtaMinCut(0.0),
	fPtMaxCut(1e8),
	fPtMinCut(0.15),
	fRunNumber(0)

{
	//
	// copy constructor - untested
	//
	((AliSpherocityUtils &) c).Copy(*this);
}//_____________________________________________________________________________
AliSpherocityUtils &AliSpherocityUtils::operator=(const AliSpherocityUtils &c)
{
	//
	// Assignment operator - untested
	//

	if (this != &c) ((AliSpherocityUtils &) c).Copy(*this);
	return *this;
}
//__________________________________________________________________________
void AliSpherocityUtils::Init(){

	fNrec = 0;

	// Define track cuts, TPC-only + TPC refit
	cout<<"-----------------------------------------------------------------------------------"<<endl;
        cout<<"-----------------------------------------------------------------------------------"<<endl;
	cout<<"INITIALIZATION OF ALI_SPHEROCITY_UTILS"<<endl;
	if(!fTrackFilterGlobal){
		cout<<"ESD analysis uses: TPC-only + TPC-refit track cuts"<<endl;	
		fTrackFilterGlobal = new AliAnalysisFilter("trackFilterTPC");
		AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts();
		esdTrackCuts->SetMinNClustersTPC(50);
		esdTrackCuts->SetMaxChi2PerClusterTPC(4);
		esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
		esdTrackCuts->SetRequireTPCRefit(kTRUE);
		esdTrackCuts->SetMaxDCAToVertexZ(3.2);
		esdTrackCuts->SetMaxDCAToVertexXY(2.4);
		esdTrackCuts->SetDCAToVertex2D(kTRUE);
		fTrackFilterGlobal->AddCuts(esdTrackCuts);
	}
	if(!fAODFilterGlobal){
		fAODFilterGlobal = 1;
		cout<<"AOD analysis uses fiter 1: TPC-only trackcuts (warning: TPC refit is not included!)"<<endl;
	}
	cout<<"Min number of particles:"<<fMinMult<<endl;
	if(fIsAbsEta)
		cout<<"Kinematic cuts:   "<<fPtMinCut<<"<pT<"<<fPtMaxCut<<" GeV/c,   "<<fEtaMinCut<<"<|eta|<"<<fEtaMaxCut<<endl;
	else
		cout<<"Kinematic cuts:   "<<fPtMinCut<<"<pT<"<<fPtMaxCut<<" GeV/c,   "<<fEtaMinCut<<"<eta<"<<fEtaMinCut<<endl;
	cout<<"Step size for spherocity calculation:"<<fSizeStep*2*TMath::Pi()/360.0<<"  radians"<<endl;
	cout<<"-----------------------------------------------------------------------------------"<<endl;
	cout<<"-----------------------------------------------------------------------------------"<<endl;

} 


//_____________________________________________________________________________
Float_t AliSpherocityUtils::MinVal( Float_t A, Float_t B ) {
	if( A < B ) {
		return A;
	}
	else {
		return B;
	}
}
//______________________________________________________________________
Float_t AliSpherocityUtils::GetEventShape( AliVEvent *event, TH1D *hphi, TH1D *heta )
{

	//Int_t lRequestedRunNumber = event->GetRunNumber();
	Double_t lreturnval = 1.0;

	//fNrec = 0;
	if (event->InheritsFrom("AliESDEvent"))
		fESDEvent = dynamic_cast<AliESDEvent *>(event);
	else if (event->InheritsFrom("AliAODEvent"))
		fAODEvent = dynamic_cast<AliAODEvent *>(event);

	lreturnval = GetSpherocity( hphi, heta );

	return lreturnval;

}
//______________________________________________________________________
Float_t AliSpherocityUtils::GetEventShapeTrue( AliStack *event, TH1D *hphi, TH1D *heta )
{

	Double_t lreturnval = 1.0;

	fMCStack = event;

	lreturnval = GetSpherocityMC( hphi, heta );

	return lreturnval;
}
//_____________________________________________________________________
Int_t AliSpherocityUtils::ReadMC( vector<Float_t> &ptArray,  vector<Float_t> &etaArray, vector<Float_t> &phiArray, TH1D * hphi, TH1D *heta ){

	ptArray.clear();
	etaArray.clear();
	phiArray.clear();

	Int_t nTracks = fMCStack->GetNtrack();
	Int_t nTrue = 0;

	for(Int_t iT = 0; iT < nTracks; iT++) {       
		//Cuts
		if(!(fMCStack->IsPhysicalPrimary(iT)))
			continue;

		TParticle* Track = fMCStack->Particle(iT);
		if(!Track)
			continue;

		Float_t eta = Track->Eta();
		Float_t pt  = Track->Pt();
		Float_t phi  = Track->Phi();

		TParticlePDG* pdgPart = Track->GetPDG();
		Double_t chargeMC = pdgPart->Charge();

		if( TMath::Abs(chargeMC) < 0.1 )
			continue;


		if(fIsAbsEta){       //cuts in pseudorapidity
			if( TMath::Abs( eta ) > fEtaMaxCut || TMath::Abs( eta ) < fEtaMinCut )
				continue;
		}
		else{
			if( eta > fEtaMaxCut || eta < fEtaMinCut )  
				continue;
		}
		//cuts in pt
		if( pt > fPtMaxCut || pt <  fPtMinCut )
			continue;

		if(hphi)
			hphi->Fill(phi);
		if(heta)
			heta->Fill(eta);

		ptArray.push_back(pt);
		etaArray.push_back(eta);
		phiArray.push_back(phi);


		nTrue++;

	} //close first loop on nTracks

	return nTrue;

}
//_____________________________________________________________________
Int_t AliSpherocityUtils::ReadAODEvent( vector<Float_t> &ptArray,  vector<Float_t> &etaArray, vector<Float_t> &phiArray, TH1D * hphi, TH1D *heta ){
	ptArray.clear();
	etaArray.clear();
	phiArray.clear();

	Int_t nTracks = fAODEvent->GetNumberOfTracks();
	Int_t nRec    = 0;

	for(Int_t iT = 0; iT < nTracks; iT++) {

		AliVTrack   *trackTmp = (AliVTrack *)fAODEvent->GetTrack(iT);
		AliAODTrack * Track  = dynamic_cast<AliAODTrack *>(trackTmp);
		if(!Track)
			continue;

		Float_t eta = Track->Eta();
		Float_t pt  = Track->Pt();
		Float_t phi  = Track->Phi();


		if(fIsAbsEta){  //cuts in pseudorapidity
			if( TMath::Abs(eta) > fEtaMaxCut || TMath::Abs(eta) < fEtaMinCut )          
				continue;
		}
		else{

			if( eta > fEtaMaxCut || eta < fEtaMinCut )
				continue;
		}
		//cuts in pt
		if( pt > fPtMaxCut || pt <  fPtMinCut )
			continue;

		//quality cuts
		//standard 2010 hard DCA
		if(!Track->TestFilterBit(fAODFilterGlobal))
			continue;


		if(hphi)
			hphi->Fill(phi);
		if(heta)
			heta->Fill(eta);

		ptArray.push_back(pt);
		etaArray.push_back(eta);
		phiArray.push_back(phi);

		nRec++;

	} //close first loop on nTracks

	return nRec;

}
//_____________________________________________________________________
Int_t AliSpherocityUtils::ReadESDEvent( vector<Float_t> &ptArray,  vector<Float_t> &etaArray, vector<Float_t> &phiArray, TH1D * hphi, TH1D *heta ){

	ptArray.clear();
	etaArray.clear();
	phiArray.clear();

	// Here I define the track cuts


	Int_t nTracks = fESDEvent->GetNumberOfTracks();
	Int_t nRec    = 0;

	for(Int_t iT = 0; iT < nTracks; iT++) {
		AliESDtrack* Track = 0;
		Track = fESDEvent->GetTrack(iT);
		if(!Track)
			continue;

		Float_t eta  = Track->Eta();
		Float_t pt   = Track->Pt();
		Float_t phi  = Track->Phi();


		if(fIsAbsEta){  //cuts in pseudorapidity
			if( TMath::Abs(eta) > fEtaMaxCut || TMath::Abs(eta) < fEtaMinCut )          
				continue;
		}
		else{

			if( eta > fEtaMaxCut || eta < fEtaMinCut )
				continue;
		}
		//cuts in pt
		if( pt > fPtMaxCut || pt <  fPtMinCut )
			continue;

		//quality cuts
		if(!fTrackFilterGlobal->IsSelected(Track))
			continue;

		if(hphi)
			hphi->Fill(phi);
		if(heta)
			heta->Fill(eta);

		ptArray.push_back(pt);
		etaArray.push_back(eta);
		phiArray.push_back(phi);

		nRec++;

	} //close first loop on nTracks

	return nRec;

}
//_____________________________________________________________________
Float_t AliSpherocityUtils::AnalyseGetSpherocity( const vector<Float_t> &pt, const vector<Float_t> &eta, const vector<Float_t> &phi ){


	Float_t spherocity = -10.0;
	Float_t pFull = 0;
	Float_t Spherocity = 2;

	//computing total pt
	Float_t sumapt = 0;
	for(Int_t i1 = 0; i1 < fNrec; ++i1){
		sumapt += pt[i1];

	}

	//Getting thrust
	for(Int_t i = 0; i < 360/(fSizeStep); ++i){
		Float_t numerador = 0;
		Float_t phiparam  = 0;
		Float_t nx = 0;
		Float_t ny = 0;
		phiparam=( (TMath::Pi()) * i * fSizeStep ) / 180; // parametrization of the angle
		nx = TMath::Cos(phiparam);            // x component of an unitary vector n
		ny = TMath::Sin(phiparam);            // y component of an unitary vector n
		for(Int_t i1 = 0; i1 < fNrec; ++i1){

			Float_t pxA = pt[i1] * TMath::Cos( phi[i1] );
			Float_t pyA = pt[i1] * TMath::Sin( phi[i1] );

			numerador += TMath::Abs( ny * pxA - nx * pyA );//product between p  proyection in XY plane and the unitary vector
		}
		pFull=TMath::Power( (numerador / sumapt),2 );
		if(pFull < Spherocity)//maximization of pFull
		{
			Spherocity = pFull;
		}
	}

	spherocity=((Spherocity)*TMath::Pi()*TMath::Pi())/4.0;


	return spherocity;

}
//_____________________________________________________________________
Float_t AliSpherocityUtils::GetSpherocity( TH1D * hphi, TH1D *heta )
{

	vector<Float_t> pt;
	vector<Float_t> eta;
	vector<Float_t> phi;

	if(fESDEvent)
		fNrec = ReadESDEvent( pt, eta, phi, hphi, heta );
	else if(fAODEvent)
		fNrec = ReadAODEvent( pt, eta, phi, hphi, heta );

	if( fNrec < fMinMult )
		return -0.5;

	Float_t spherocity = AnalyseGetSpherocity( pt, eta, phi ); 


	return spherocity;

}
//_____________________________________________________________________
Float_t AliSpherocityUtils::GetSpherocityMC( TH1D * hphi, TH1D *heta )
{
	vector<Float_t> pt;
	vector<Float_t> eta;
	vector<Float_t> phi;

	fNrec = ReadMC(pt, eta, phi, hphi, heta); 
	if( fNrec < fMinMult )
		return -0.5;

	Float_t spherocity = AnalyseGetSpherocity( pt, eta, phi ); 

	return spherocity;
}

