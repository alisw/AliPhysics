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
   Class AliTransverseEventShape  
   It allow to calculate Sphericity and Spherocity variables

 */
//           Please report any Bugs to:
//*************************************************************************
//   author: Antonio Ortiz Velasquez  ( antonio.ortiz@nucleares.unam.mx )
//           Eleazar Cuautle Flores   ( ecuautle@nucleares.unam.mx )
/**************************************************************************


  First version, June 4 2015
  02/07/2015, Antonio: adding MC true event shape calculation
  16/07/2015, Antonio: AODs are supported

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
#include "AliTransverseEventShape.h"
#include <TFile.h>
#include "AliAODHeader.h"
// STL includes
#include <iostream>
using namespace std;


ClassImp(AliTransverseEventShape)

	//______________________________________________________________________
AliTransverseEventShape::AliTransverseEventShape():TObject(),
	fAnalysisMC(kFALSE),
	fESDEvent(0),
	fAODEvent(0),
	fMCStack(0),
	fNrec(0),
	fUseHybrid(0),
	fTrackFilterHybrid1(0),
	fTrackFilterHybrid2(0),
	fTrackFilterGlobal(0),
        fAODFilterHybrid1(0),
        fAODFilterHybrid2(0),
        fAODFilterGlobal(0),
	fMinMultESA(0),
	fSizeStepESA(0),
	fIsAbsEtaESA(0),
	fEtaMaxCutESA(0),
	fEtaMinCutESA(0),
	fPtMaxCutESA(0),
	fPtMinCutESA(0),
	fRunNumber(0),
	fhetaSo(0),
	fhphiSo(0),
	fhptSo(0),
	fhetaSt(0),
	fhphiSt(0),
	fhptSt(0),
	fhetaSoMC(0),
	fhphiSoMC(0),
	fhptSoMC(0),
	fhetaStMC(0),
	fhphiStMC(0),
	fhptStMC(0)

{
	// Default contructor
}


//_____________________________________________________________________________
AliTransverseEventShape::AliTransverseEventShape(const AliTransverseEventShape &c) : 
	TObject(c),
	fESDEvent(0),
	fAODEvent(0),
	fMCStack(0),
	fNrec(0),
	fUseHybrid(0),
	fAnalysisMC(kFALSE),
	fTrackFilterHybrid1(0),
	fTrackFilterHybrid2(0),
	fTrackFilterGlobal(0),
        fAODFilterHybrid1(0),
        fAODFilterHybrid2(0),
        fAODFilterGlobal(0),
	fMinMultESA(0),
	fSizeStepESA(0),
	fIsAbsEtaESA(0),
	fEtaMaxCutESA(0),
	fEtaMinCutESA(0),
	fPtMaxCutESA(0),
	fPtMinCutESA(0),
	fRunNumber(0),
	fhetaSo(0),
	fhphiSo(0),
	fhptSo(0),
	fhetaSt(0),
	fhphiSt(0),
	fhptSt(0),
	fhetaSoMC(0),
	fhphiSoMC(0),
	fhptSoMC(0),
	fhetaStMC(0),
	fhphiStMC(0),
	fhptStMC(0)

{
	//
	// copy constructor - untested
	//
	((AliTransverseEventShape &) c).Copy(*this);
}//_____________________________________________________________________________
AliTransverseEventShape &AliTransverseEventShape::operator=(const AliTransverseEventShape &c)
{
	//
	// Assignment operator - untested
	//

	if (this != &c) ((AliTransverseEventShape &) c).Copy(*this);
	return *this;
}
//__________________________________________________________________________
void AliTransverseEventShape::Init(){

	fNrec = 0;

	fhetaSo = new TH1D("fhetaSo","ESD;#eta; entries",30,-1.5,1.5);
	fhphiSo = new TH1D("fhphiSo","ESD;#phi (rad); entries", 64, 0, 2*TMath::Pi());
	fhptSo  = new TH1D("fhptSo","ESD;#it{p}_{T} (GeV/#it{c}); entries", 1000, 0, 100);

	fhetaSt = new TH1D("fhetaSt","ESD;#eta; entries",30,-1.5,1.5);
	fhphiSt = new TH1D("fhphiSt","ESD;#phi (rad); entries", 64, 0, 2*TMath::Pi());
	fhptSt  = new TH1D("fhptSt","ESD;#it{p}_{T} (GeV/#it{c}); entries", 1000, 0, 100);

	if( fAnalysisMC ){
		fhetaSoMC = new TH1D("fhetaSoMC","MC;#eta; entries",30,-1.5,1.5);
		fhphiSoMC = new TH1D("fhphiSoMC","MC;#phi (rad); entries", 64, 0, 2*TMath::Pi());
		fhptSoMC  = new TH1D("fhptSoMC","MC;#it{p}_{T} (GeV/#it{c}); entries", 1000, 0, 100);

		fhetaStMC = new TH1D("fhetaStMC","MC;#eta; entries",30,-1.5,1.5);
		fhphiStMC = new TH1D("fhphiStMC","MC;#phi (rad); entries", 64, 0, 2*TMath::Pi());
		fhptStMC  = new TH1D("fhptStMC","MC;#it{p}_{T} (GeV/#it{c}); entries", 1000, 0, 100);
	}
} 


//_____________________________________________________________________________
Float_t AliTransverseEventShape::MinVal( Float_t A, Float_t B ) {
	if( A < B ) {
		return A;
	}
	else {
		return B;
	}
}
//______________________________________________________________________
Float_t AliTransverseEventShape::GetEventShape( AliVEvent *event, TString lMethod, Bool_t fillHist )
{

	Int_t lRequestedRunNumber = event->GetRunNumber();
	Double_t lreturnval = 1.0;

	//fNrec = 0;
	if (event->InheritsFrom("AliESDEvent"))
		fESDEvent = dynamic_cast<AliESDEvent *>(event);
	else if (event->InheritsFrom("AliAODEvent"))
		fAODEvent = dynamic_cast<AliAODEvent *>(event);

	if ( lMethod == "SO" ) lreturnval = GetSpherocity( fillHist );
	if ( lMethod == "ST" ) lreturnval = GetSphericity( fillHist );

	return lreturnval;

}
//______________________________________________________________________
Float_t AliTransverseEventShape::GetEventShapeTrue( AliStack *event, TString lMethod, Bool_t fillHist )
{

	Double_t lreturnval = 1.0;

	fMCStack = event;

	if ( lMethod == "SO" ) lreturnval = GetSpherocityMC( fillHist );
	if ( lMethod == "ST" ) lreturnval = GetSphericityMC( fillHist );

	return lreturnval;
}
//_____________________________________________________________________
Int_t AliTransverseEventShape::ReadMC( vector<Float_t> &ptArray,  vector<Float_t> &etaArray, vector<Float_t> &phiArray ){

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


		if(fIsAbsEtaESA){       //cuts in pseudorapidity
			if( TMath::Abs( eta ) > fEtaMaxCutESA || TMath::Abs( eta ) < fEtaMinCutESA )
				continue;
		}
		else{
			if( eta > fEtaMaxCutESA || eta < fEtaMinCutESA )  
				continue;
		}
		//cuts in pt
		if( pt > fPtMaxCutESA || pt <  fPtMinCutESA )
			continue;

		ptArray.push_back(pt);
		etaArray.push_back(eta);
		phiArray.push_back(phi);

		nTrue++;

	} //close first loop on nTracks

	return nTrue;

}
//_____________________________________________________________________
Int_t AliTransverseEventShape::ReadAODEvent( vector<Float_t> &ptArray,  vector<Float_t> &etaArray, vector<Float_t> &phiArray ){

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


		if(fIsAbsEtaESA){  //cuts in pseudorapidity
			if( TMath::Abs(eta) > fEtaMaxCutESA || TMath::Abs(eta) < fEtaMinCutESA )          
				continue;
		}
		else{

			if( eta > fEtaMaxCutESA || eta < fEtaMinCutESA )
				continue;
		}
		//cuts in pt
		if( pt > fPtMaxCutESA || pt <  fPtMinCutESA )
			continue;

		//quality cuts
		//following this presentation by Martha: https://twiki.cern.ch/twiki/pub/ALICE/PWGPPAODTrackCuts/FilterBitProposal.pdf
		if(fUseHybrid){

			//Bool_t ishy1 = Track->TestFilterBit(256);
			//Bool_t ishy2 = Track->TestFilterBit(512);
                        Bool_t ishy1 = Track->TestFilterBit(fAODFilterHybrid1);
                        Bool_t ishy2 = Track->TestFilterBit(fAODFilterHybrid2);

			if( !(ishy1 || ishy2) )
				continue;

		}else{
			//standard 2010 hard DCA
			if(!Track->TestFilterBit(fAODFilterGlobal))
				continue;
		}

		ptArray.push_back(pt);
		etaArray.push_back(eta);
		phiArray.push_back(phi);

		nRec++;

	} //close first loop on nTracks

	return nRec;

}
//_____________________________________________________________________
Int_t AliTransverseEventShape::ReadESDEvent( vector<Float_t> &ptArray,  vector<Float_t> &etaArray, vector<Float_t> &phiArray ){

	ptArray.clear();
	etaArray.clear();
	phiArray.clear();

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


		if(fIsAbsEtaESA){  //cuts in pseudorapidity
			if( TMath::Abs(eta) > fEtaMaxCutESA || TMath::Abs(eta) < fEtaMinCutESA )          
				continue;
		}
		else{

			if( eta > fEtaMaxCutESA || eta < fEtaMinCutESA )
				continue;
		}
		//cuts in pt
		if( pt > fPtMaxCutESA || pt <  fPtMinCutESA )
			continue;

		//quality cuts
		if(!fUseHybrid){//golden track cuts
			if(!fTrackFilterGlobal->IsSelected(Track))
				continue;
		}
		else{//hybrid track cuts
			Bool_t cutset1 = kFALSE;
			Bool_t cutset2 = kFALSE;
			cutset1 = fTrackFilterHybrid1->IsSelected(Track);
			cutset2 = (!fTrackFilterGlobal->IsSelected(Track)) && fTrackFilterHybrid2->IsSelected(Track);
			if(!(cutset1 || cutset2))
				continue;
		}

		ptArray.push_back(pt);
		etaArray.push_back(eta);
		phiArray.push_back(phi);

		nRec++;

	} //close first loop on nTracks

	return nRec;

}
//_____________________________________________________________________
Float_t AliTransverseEventShape::AnalyseGetSphericity( Bool_t fillHist, const vector<Float_t> &pt, const vector<Float_t> &eta, const vector<Float_t> &phi ){


	Float_t sphericity = -10;
	Float_t s00=0;
	Float_t s01=0;
	Float_t s11=0;
	Float_t totalpt=0;

	for(Int_t i1 = 0; i1 < fNrec; ++i1){


		//Fill QA histos
		if(fillHist){
			fhetaSt->Fill(eta[i1]);
			fhphiSt->Fill(phi[i1]);
			fhptSt->Fill(pt[i1]);
		}
		Float_t px = pt[i1] * TMath::Cos( phi[i1] );
		Float_t py = pt[i1] * TMath::Sin( phi[i1] );

		s00 += (px * px)/pt[i1];
		s01 += (py * px)/pt[i1];
		s11 += (py * py)/pt[i1];
		totalpt += pt[i1];
	}


	Double_t S00=s00/totalpt;
	Double_t S01=s01/totalpt;
	Double_t S11=s11/totalpt;

	Float_t lambda1=((S00+S11)+TMath::Sqrt((S00+S11)*(S00+S11)-4*(S00*S11-S01*S01)))/2;
	Float_t lambda2=((S00+S11)-TMath::Sqrt((S00+S11)*(S00+S11)-4*(S00*S11-S01*S01)))/2;
	if((lambda2==0)&&(lambda1==0))
		sphericity=0;
	if(lambda1+lambda2!=0)
		sphericity=2*TMath::Min( lambda1,lambda2 )/( lambda1+lambda2 );

	return sphericity;

}


//_____________________________________________________________________
Float_t AliTransverseEventShape::AnalyseGetSpherocity( Bool_t fillHist, const vector<Float_t> &pt, const vector<Float_t> &eta, const vector<Float_t> &phi ){


	Float_t spherocity = -10.0;
	Float_t pFull = 0;
	Float_t Spherocity = 2;

	//computing total pt
	Float_t sumapt = 0;
	for(Int_t i1 = 0; i1 < fNrec; ++i1){
		sumapt += pt[i1];

		//Fill QA histos
		if(fillHist){
			fhetaSo->Fill(eta[i1]);
			fhphiSo->Fill(phi[i1]);
			fhptSo->Fill(pt[i1]);
		}

	}

	//Getting thrust
	for(Int_t i = 0; i < 360/(fSizeStepESA); ++i){
		Float_t numerador = 0;
		Float_t phiparam  = 0;
		Float_t nx = 0;
		Float_t ny = 0;
		phiparam=( (TMath::Pi()) * i * fSizeStepESA ) / 180; // parametrization of the angle
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
Float_t AliTransverseEventShape::GetSpherocity( Bool_t fillHist )
{

	vector<Float_t> pt;
	vector<Float_t> eta;
	vector<Float_t> phi;

	if(fESDEvent)
		fNrec = ReadESDEvent(pt, eta, phi);
	else if(fAODEvent)
		fNrec = ReadAODEvent(pt, eta, phi);

	if( fNrec < fMinMultESA )
		return -0.5;

	Float_t spherocity = AnalyseGetSpherocity( fillHist, pt, eta, phi ); 

	return spherocity;

}
//____________________________________________________________________
Float_t AliTransverseEventShape::GetSphericity( Bool_t fillHist  )
{


	vector<Float_t> pt;
	vector<Float_t> eta;
	vector<Float_t> phi;

	if(fESDEvent)
		fNrec = ReadESDEvent(pt, eta, phi);
	else if(fAODEvent)
		fNrec = ReadAODEvent(pt, eta, phi);

	if( fNrec < fMinMultESA )
		return -0.5;

	Float_t sphericity = AnalyseGetSphericity( fillHist, pt, eta, phi );

	return sphericity;


}
//____________________________________________________________________
Float_t AliTransverseEventShape::GetSphericityMC( Bool_t fillHist )
{

	vector<Float_t> pt;
	vector<Float_t> eta;
	vector<Float_t> phi;

	fNrec = ReadMC(pt, eta, phi);
	if( fNrec < fMinMultESA )
		return -0.5;

	Float_t sphericity = AnalyseGetSphericity( fillHist, pt, eta, phi );

	return sphericity;

}
//_____________________________________________________________________
Float_t AliTransverseEventShape::GetSpherocityMC( Bool_t fillHist )
{
	vector<Float_t> pt;
	vector<Float_t> eta;
	vector<Float_t> phi;

	fNrec = ReadMC(pt, eta, phi); 
	if( fNrec < fMinMultESA )
		return -0.5;

	Float_t spherocity = AnalyseGetSpherocity( fillHist, pt, eta, phi ); 

	return spherocity;
}
//________________________________________________________
TH1D * AliTransverseEventShape::GetHistData( Int_t bin_histo, TString lMethod ){

	if ( lMethod == "SO" ){
		switch(bin_histo){
			case 0:
				return fhetaSo;
				break;
			case 1:
				return fhphiSo;
				break;
			case 2:
				return fhptSo;
				break;
			default:
				return 0;
		}
	}
	else if ( lMethod == "ST" ){
		switch(bin_histo){
			case 0:
				return fhetaSt;
				break;
			case 1:
				return fhphiSt;
				break;
			case 2:
				return fhptSt;
				break;
			default:
				return 0;
		}
	}



	return 0;
}
//__________________________________________________________
void AliTransverseEventShape::SaveHistosSo( const char* folder ){

	if(!fhetaSo)
		return;

	if (folder)
	{
		gDirectory->mkdir(folder);
		gDirectory->cd(folder);
	}

	fhetaSo->Write();
	fhphiSo->Write();
	fhptSo->Write();
	if(fAnalysisMC){

		fhetaSoMC->Write();
		fhphiSoMC->Write();
		fhptSoMC->Write();

	}
	if (folder)
		gDirectory->cd("..");

}
//__________________________________________________________
void AliTransverseEventShape::SaveHistosSt( const char* folder ){

	if(!fhetaSt)
		return;

	if (folder)
	{
		gDirectory->mkdir(folder);
		gDirectory->cd(folder);
	}

	fhetaSt->Write();
	fhphiSt->Write();
	fhptSt->Write();
	if(fAnalysisMC){

		fhetaStMC->Write();
		fhphiStMC->Write();
		fhptStMC->Write();

	}
	if (folder)
		gDirectory->cd("..");

}

